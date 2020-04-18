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
from itertools import product, repeat
import os

APS = [
    CircularAperture,
    CircularAnnulus,
    EllipticalAperture,
    EllipticalAnnulus,
    RectangularAperture,
    RectangularAnnulus,
]
params = [(3,), (3, 5), (3, 3, 0), (3, 5, 4, 0), (3, 5, 0), (3, 5, 4, 0)]

Ns = [1, 10, 50, 100]
methods = ["exact", "center", "subpixel"]

data = np.random.randn(512, 512) + 10

inputs = list(product(Ns, methods, zip(APS, params)))

rows = []
for input_ in tqdm(inputs):
    N, method, (A, par) = input_
    if method == "exact" and (A == RectangularAnnulus or A == RectangularAperture):
        continue
    aps = A(list(repeat((255, 255), N)), *par)
    m = method if method != "subpixel" else "subpixel-5"
    kwargs = {} if method != "subpixel" else {"subpixels": 5}
    t0 = datetime.now()
    t = aperture_photometry(data, aps, method=method, **kwargs)
    t1 = datetime.now()
    time = (t1 - t0).total_seconds()
    rows.append((N, A.__name__, m, time))

df = pd.DataFrame(rows, columns=["N", "aperture", "method", "time"])

path = os.path.dirname(__file__)
df.to_csv(os.path.join(path, "python_apertures.csv"), index=False)
