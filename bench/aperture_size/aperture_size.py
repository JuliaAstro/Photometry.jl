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
from tqdm import tqdm, trange
from itertools import repeat
import os

rng = np.random.seed(11256)

data = np.random.randn(512, 512) + 10

rows = []
for r in trange(1, 201, 5):
    ap = CircularAperture((255.5, 255.5), r)

    ts = []
    for i in range(5):
        t0 = datetime.now()
        t = aperture_photometry(data, ap, method="exact")
        t1 = datetime.now()
        time = (t1 - t0).total_seconds()
        ts.append(time)

    time = sum(ts) / 5
    rows.append((r, time))

df = pd.DataFrame(rows, columns=["r", "time"])

path = os.path.dirname(__file__)
df.to_csv(os.path.join(path, "python_aperture_size.csv"), index=False)

## Ellipse

rows = []
for r in trange(1, 201, 5):
    ap = EllipticalAperture((255.5, 255.5), r, r, 20)

    ts = []
    for i in range(5):
        t0 = datetime.now()
        t = aperture_photometry(data, ap, method="exact")
        t1 = datetime.now()
        time = (t1 - t0).total_seconds()
        ts.append(time)

    time = sum(ts) / 5
    rows.append((r, time))

df = pd.DataFrame(rows, columns=["r", "time"])

path = os.path.dirname(__file__)
df.to_csv(os.path.join(path, "python_aperture_size-ellipse.csv"), index=False)
