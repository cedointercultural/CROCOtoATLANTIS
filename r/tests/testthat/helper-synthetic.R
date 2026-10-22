# Synthetic BGM and CROCO file: two boxes of 1 x 2 km side by side, 100 m deep,
# and a CROCO grid of 6 x 6 cells of 500 m (4 x 4 interior cells, closed ghost
# ring) with a uniform current of 1 m/s from box 0 to box 1.

skip_if_no_python <- function() {
  skip_if_not(croco2atlantis_available(), "the Python package croco2atlantis is not available")
}

BOTTOMS <- c(18, 54, 270)
PROJ <- "+proj=laea +lat_0=26 +lon_0=-110 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"

synthetic_bgm <- function(dir) {
  box <- function(i, x0, other) {
    c(sprintf("box%d.label\tBox%d", i, i),
      sprintf("box%d.inside\t%d 1000", i, x0 + 500),
      sprintf("box%d.nconn\t1", i),
      sprintf("box%d.iface\t0", i),
      sprintf("box%d.ibox\t%d", i, other),
      sprintf("box%d.botz\t-100.0", i),
      sprintf("box%d.area\t2000000", i),
      sprintf("box%d.vertmix\tNaN", i),
      sprintf("box%d.horizmix\tNaN", i),
      sprintf("box%d.vert %d %d", i, c(x0, x0 + 1000, x0 + 1000, x0, x0), c(0, 0, 2000, 2000, 0)))
  }
  lines <- c("# test", paste("projection", PROJ), "nbox 2", "nface 1", "maxwcbotz -100.0",
             box(0, 0, 1), box(1, 1000, 0),
             "face0.p1 1000 0", "face0.p2 1000 2000", "face0.length 2000", "face0.cs 0 1", "face0.lr 0 1")
  path <- file.path(dir, "test.bgm")
  writeLines(lines, path)
  path
}

synthetic_croco <- function(dir) {
  path <- file.path(dir, "croco_avg.nc")
  reticulate::py_run_string(sprintf('
import numpy as np, xarray as xr
from pyproj import Transformer
n, d, ns, nt = 6, 500.0, 4, 2
c = np.arange(n) * d - d / 2  # centres from -250 m; interior cells 250..1750 m
y, x = np.meshgrid(c, c, indexing="ij")
lon, lat = Transformer.from_crs("%s", "EPSG:4326", always_xy=True).transform(x, y)
mask = np.zeros((n, n)); mask[1:-1, 1:-1] = 1
s_w = np.linspace(-1, 0, ns + 1); s_rho = (s_w[:-1] + s_w[1:]) / 2
temp = np.where(x < 1000, 15.0, 20.0)[None, None] * np.ones((nt, ns, n, n))
rho = ("time", "s_rho", "eta_rho", "xi_rho")
two = ("eta_rho", "xi_rho")
xr.Dataset({
    "lon_rho": (two, lon), "lat_rho": (two, lat), "h": (two, np.full((n, n), 100.0)),
    "mask_rho": (two, mask), "angle": (two, np.zeros((n, n))),
    "pm": (two, np.full((n, n), 1 / d)), "pn": (two, np.full((n, n), 1 / d)),
    "Vtransform": ((), 2.0), "hc": ((), 10.0),
    "s_rho": (("s_rho",), s_rho), "s_w": (("s_w",), s_w),
    "Cs_rho": (("s_rho",), s_rho), "Cs_w": (("s_w",), s_w),
    "time": (("time",), np.array([0.0, 43200.0])),
    "zeta": (("time",) + two, np.zeros((nt, n, n))),
    "temp": (rho, temp), "salt": (rho, np.full((nt, ns, n, n), 35.0)),
    "u": (("time", "s_rho", "eta_rho", "xi_u"), np.ones((nt, ns, n, n - 1))),
    "v": (("time", "s_rho", "eta_v", "xi_rho"), np.zeros((nt, ns, n - 1, n))),
}).to_netcdf(r"%s")
', PROJ, path))
  path
}
