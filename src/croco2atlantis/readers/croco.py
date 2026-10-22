"""Read CROCO output files.

The conversion uses the average files (`avg`): in a run with tides, instantaneous
history files (`his`) every 12 h alias the semidiurnal tide into slow, spurious
signals. Both kinds have the same layout, so everything here works on either.

Works on the raw NetCDF files (dimensions `time`, `s_rho`, `eta_rho`, `xi_rho`, ...)
without `croco_pytools`. Ghost points are kept: the rho grid stays (eta_rho, xi_rho)
and the land/ghost cells are removed with the mask. Land cells hold 0 in CROCO, not
NaN, so every average must be masked.

Level order: index 0 is the bottom and the last index the surface.
"""

from dataclasses import dataclass

import numpy as np
import pandas as pd
import xarray as xr


@dataclass
class CrocoGrid:
    """Horizontal grid on rho points, all arrays shaped (eta_rho, xi_rho)."""

    lon: np.ndarray
    lat: np.ndarray
    h: np.ndarray  # bathymetry, m, positive
    mask: np.ndarray  # bool, True on water, False on the ghost ring
    angle: np.ndarray  # angle between the xi axis and east, rad
    pm: np.ndarray  # 1 / cell size in xi, 1/m
    pn: np.ndarray  # 1 / cell size in eta, 1/m
    open_boundary: np.ndarray | None = None  # bool, water cells of the ghost ring

    @property
    def cell_area(self):
        """Area of each cell in m2."""
        return 1.0 / (self.pm * self.pn)


@dataclass
class VerticalParams:
    """Parameters of the CROCO s-coordinate."""

    vtransform: int
    hc: float
    s_rho: np.ndarray
    s_w: np.ndarray
    cs_rho: np.ndarray
    cs_w: np.ndarray


def open_output(path):
    """Open one CROCO average (or history) file lazily. Call once per file."""
    return xr.open_dataset(path, decode_times=False)


def read_grid(ds):
    """Horizontal grid of an open history file.

    The outer ring of rho points are ghost points of the model; they are marked as
    land in `mask` so they never enter an average (13909 water cells remain in the
    Gulf of California grid, the same as `xcroco`). The ghost cells that are water
    are the open boundaries of the model and are kept in `open_boundary`, for the
    transport through them.
    """
    water = ds["mask_rho"].values == 1
    mask = water.copy()
    mask[[0, -1], :] = False
    mask[:, [0, -1]] = False
    return CrocoGrid(
        lon=ds["lon_rho"].values,
        lat=ds["lat_rho"].values,
        h=ds["h"].values.astype(float),
        mask=mask,
        angle=ds["angle"].values.astype(float),
        pm=ds["pm"].values.astype(float),
        pn=ds["pn"].values.astype(float),
        open_boundary=water & ~mask,
    )


def read_vertical_params(ds):
    """s-coordinate parameters of an open history file."""
    def first(name):
        v = np.asarray(ds[name].values, dtype=float)
        return v.reshape(-1, v.shape[-1])[0] if v.ndim > 1 else v

    vtransform = int(np.asarray(ds["Vtransform"].values).ravel()[0])
    if vtransform not in (1, 2):
        raise ValueError(f"Vtransform must be 1 or 2, found {vtransform}")
    return VerticalParams(
        vtransform=vtransform,
        hc=float(np.asarray(ds["hc"].values).ravel()[0]),
        s_rho=first("s_rho"),
        s_w=first("s_w"),
        cs_rho=first("Cs_rho"),
        cs_w=first("Cs_w"),
    )


def read_times(ds, origin):
    """Dates of the records. `origin` is the date the `time` seconds count from.

    The files only say "time since initialization", so the origin is not in them and
    must be given (for the Gulf of California run it is consistent with 2000-01-01).
    In average files the time of a record is the centre of its averaging window.
    """
    seconds = ds["time"].values
    return pd.Timestamp(origin) + pd.to_timedelta(seconds, unit="s")


def s_depths(h, zeta, params, where="w"):
    """Depth of the s-levels, in metres, negative below the surface.

    `h` and `zeta` are arrays of the same shape (the horizontal dimensions). Returns
    an array shaped (n_levels,) + h.shape, with index 0 at the bottom. `where` is
    "w" for the interfaces (n_s + 1 levels) or "rho" for the level centres.
    """
    if where not in ("w", "rho"):
        raise ValueError(f'where must be "w" or "rho", found {where!r}')
    s, cs = (params.s_w, params.cs_w) if where == "w" else (params.s_rho, params.cs_rho)
    pad = (slice(None),) + (None,) * np.ndim(h)
    s, cs = s[pad], cs[pad]
    h = np.asarray(h, dtype=float)
    zeta = np.asarray(zeta, dtype=float)
    hc = params.hc
    if params.vtransform == 2:
        z0 = (hc * s + h * cs) / (hc + h)
        return zeta + (zeta + h) * z0
    z0 = hc * s + (h - hc) * cs
    return z0 + zeta * (1.0 + z0 / h)
