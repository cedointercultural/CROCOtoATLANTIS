# CROCO output.

#' Open a CROCO output file
#'
#' Opens one average (`avg`) or history (`his`) file lazily with xarray. The
#' conversion should use the average files: in a run with tides, 12 h snapshots
#' alias the semidiurnal tide.
#'
#' @param path Path of the NetCDF file.
#' @return The Python xarray dataset.
#' @export
croco_open <- function(path) {
  py_mod("readers.croco")$open_output(normalizePath(path, mustWork = TRUE))
}

#' Horizontal grid of a CROCO file
#'
#' The ghost ring of the model is masked out; its water cells are the open
#' boundaries.
#'
#' @param ds A dataset from [croco_open()].
#' @return The Python `CrocoGrid` object.
#' @export
croco_grid <- function(ds) py_mod("readers.croco")$read_grid(ds)

#' Parameters of the s-coordinate of a CROCO file
#'
#' @param ds A dataset from [croco_open()].
#' @return The Python `VerticalParams` object.
#' @export
croco_vertical_params <- function(ds) py_mod("readers.croco")$read_vertical_params(ds)

#' Dates of the records of a CROCO file
#'
#' The files count time from the start of the run, which they do not record, so
#' the origin must be given. In average files the date of a record is the centre
#' of its averaging window.
#'
#' @param ds A dataset from [croco_open()].
#' @param origin Date the `time` seconds count from, for example `"2000-01-01"`.
#' @return A `POSIXct` vector in UTC.
#' @export
croco_times <- function(ds, origin) {
  t <- py_mod("readers.croco")$read_times(ds, as.character(origin))
  as.POSIXct(reticulate::py_to_r(t$strftime("%Y-%m-%d %H:%M:%S")$tolist()), tz = "UTC")
}

#' Depth of the s-levels at one record
#'
#' @param ds A dataset from [croco_open()].
#' @param grid A grid from [croco_grid()].
#' @param params Parameters from [croco_vertical_params()].
#' @param time Record, from 1.
#' @param where `"w"` for the level interfaces or `"rho"` for the level centres.
#' @return An array (level, eta, xi) in metres, negative below the surface, with
#'   level 1 at the bottom (the CROCO order).
#' @export
croco_s_depths <- function(ds, grid, params, time = 1, where = c("w", "rho")) {
  where <- match.arg(where)
  reticulate::py_to_r(s_depths_py(ds, grid, params, time, where))
}

s_depths_py <- function(ds, grid, params, time, where = "w") {
  zeta <- ds$zeta$isel(time = py_index(time))$values
  py_mod("readers.croco")$s_depths(grid$h, zeta, params, where)
}

croco_records <- function(ds) seq_len(reticulate::py_len(ds$time))

croco_field <- function(ds, name, time) {
  reticulate::py_get_item(ds, name)$isel(time = py_index(time))$values
}
