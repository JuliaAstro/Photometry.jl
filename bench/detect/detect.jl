using Photometry
using Chairmarks
using CSV
using DataFramesMeta
using ProgressMeter
using Random
using Statistics
using PythonCall

# ── Shared parameters ────────────────────────────────────────────────────────
const NSIGMA   = 3.0
const IMG_SIZE = 512

rng  = Random.seed!(11256)
data = randn(rng, IMG_SIZE, IMG_SIZE)

# ── Julia benchmark ───────────────────────────────────────────────────────────
julia_rows = @showprogress "Julia (PeakMesh)" map(3:2:51) do box
    alg   = PeakMesh(box_size = (box, box), nsigma = NSIGMA)
    bench = @b extract_sources($alg, $data)
    (
        box_size = box,
        nt       = Threads.nthreads(),
        median_s = bench.time,
        allocs   = bench.allocs,
    )
end

# Write per-thread file
thread_path = joinpath(@__DIR__, "julia_extract_sources_nt$(Threads.nthreads()).csv")
CSV.write(thread_path, DataFrame(julia_rows))
@info "Thread results written" path=thread_path

# Merge all per-thread files into combined CSV
julia_path  = joinpath(@__DIR__, "julia_extract_sources.csv")
thread_files = filter(
    f -> occursin(r"julia_extract_sources_nt\d+\.csv", f),
    readdir(@__DIR__, join=true)
)

if !isempty(thread_files)
    combined = mapreduce(f -> CSV.read(f, DataFrame), vcat, thread_files)
    sort!(combined, [:nt, :box_size])
    CSV.write(julia_path, combined)
    @info "Combined results written" path=julia_path nthreads=unique(combined.nt)
else
    @warn "No per-thread files found, skipping merge"
end

# ── Python benchmark ──────────────────────────────────────────────────────────
python_path = joinpath(@__DIR__, "python_extract_sources.csv")

if !isfile(python_path)
    py_script = """
import numpy as np
import pandas as pd
import os
import time
from photutils.detection import find_peaks

NSIGMA   = $(NSIGMA)
IMG_SIZE = $(IMG_SIZE)
N_TRIALS = 7
WARMUP   = 2

rng  = np.random.default_rng(11256)
data = rng.standard_normal((IMG_SIZE, IMG_SIZE))

threshold = NSIGMA * data.std()

rows = []
for box in range(3, 52, 2):
    times = []
    for i in range(WARMUP + N_TRIALS):
        t0 = time.perf_counter()
        tbl = find_peaks(data, threshold, box_size=box)
        t1 = time.perf_counter()
        if i >= WARMUP:
            times.append(t1 - t0)

    rows.append({
        "box_size" : box,
        "median_s" : float(np.median(times)),
        "min_s"    : float(np.min(times)),
        "n_sources": 0 if tbl is None else len(tbl),
    })
    print(f"box={box:3d}  median={np.median(times)*1e3:.2f} ms  "
          f"sources={rows[-1]['n_sources']}")

df = pd.DataFrame(rows)
df.to_csv("$python_path", index=False)
print("Written: $python_path")
"""

    @info "Running Python find_peaks benchmark…"
    pyexec(py_script, Main)
    @info "Python benchmark complete"
else
    @info "Skipping Python benchmark, file already exists" path=python_path
end
