from photutils import (
    CircularAperture,
    CircularAnnulus,
    EllipticalAperture,
    EllipticalAnnulus,
    RectangularAperture,
    RectangularAnnulus,
    aperture_photometry,
)
from datetime import datetime
import pandas as pd
import numpy as np
from tqdm import tqdm
from itertools import repeat
import os

rng = np.random.seed(11256)

data = np.random.randn(512, 512) + 10

rows = []
for N in tqdm([1, 10, 50, 100, 200, 400, 500, 1000, 2000]):
    ys = data.shape[0] * np.random.rand(N)
    xs = data.shape[1] * np.random.rand(N)
    pos = [(x, y) for x, y in zip(xs, ys)]
    aps = CircularAperture(pos, 10)

    ts = []
    for _ in range(5):
        t0 = datetime.now()
        t = aperture_photometry(data, aps, method="exact")
        t1 = datetime.now()
        ts.append((t1 - t0).total_seconds())
        
    time = sum(ts) / 5
    rows.append((N, time))

df = pd.DataFrame(rows, columns=["N", "time"])

path = os.path.dirname(__file__)
df.to_csv(os.path.join(path, "python_num_apertures.csv"), index=False)

# Ellipse

rows = []
for N in tqdm([1, 10, 50, 100, 200, 400, 500, 1000, 2000]):
    ys = data.shape[0] * np.random.rand(N)
    xs = data.shape[1] * np.random.rand(N)
    pos = [(x, y) for x, y in zip(xs, ys)]
    aps = EllipticalAperture(pos, 10, 10, 20)

    ts = []
    for _ in range(5):
        t0 = datetime.now()
        t = aperture_photometry(data, aps, method="exact")
        t1 = datetime.now()
        ts.append((t1 - t0).total_seconds())
        
    time = sum(ts) / 5
    rows.append((N, time))

df = pd.DataFrame(rows, columns=["N", "time"])

path = os.path.dirname(__file__)
df.to_csv(os.path.join(path, "python_num_apertures-ellipse.csv"), index=False)
