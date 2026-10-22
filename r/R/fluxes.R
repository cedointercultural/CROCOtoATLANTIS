# Water transport between boxes.

#' Cell faces between boxes
#'
#' Finds the faces of the CROCO cells (C-grid) whose two cells belong to
#' different boxes; the transport between boxes is taken on them. Call it once
#' per cell map.
#'
#' @param cellmap A cell map from [assign_cells()].
#' @param grid The grid used for `cellmap`, from [croco_grid()].
#' @return The Python `CellFaces` object.
#' @export
cell_faces <- function(cellmap, grid) py_mod("fluxes")$cell_faces(cellmap, grid)

transport_py <- function(ds, grid, params, faces, bottoms, t) {
  py_mod("fluxes")$face_transport(
    faces, croco_field(ds, "u", t), croco_field(ds, "v", t),
    s_depths_py(ds, grid, params, t), np_float(bottoms)
  )
}

face_pairs <- function(faces) {
  matrix(as.integer(reticulate::py_to_r(faces$pairs)), ncol = 2)
}

#' Water transport between boxes per layer
#'
#' Layers here are the standard ones, from the free surface to the bottom of the
#' water column, without truncating them at the depth of each box.
#'
#' @inheritParams box_layer_means
#' @param faces Cell faces from [cell_faces()].
#' @return A data frame with one row per record, pair of boxes and layer: `time`,
#'   `from`, `to`, `layer` and `transport` (m3/s, positive from `from` to `to`).
#'   Box -1 is outside the boxes (the open boundaries of CROCO).
#' @export
box_transport <- function(ds, grid, params, faces, bottoms, time = NULL) {
  if (is.null(time)) time <- croco_records(ds)
  pairs <- face_pairs(faces)
  nlayer <- length(bottoms)
  rows <- lapply(time, function(t) {
    q <- reticulate::py_to_r(transport_py(ds, grid, params, faces, bottoms, t))
    data.frame(
      time = as.integer(t),
      from = rep(pairs[, 1], nlayer),
      to = rep(pairs[, 2], nlayer),
      layer = rep(seq_len(nlayer), each = nrow(pairs)),
      transport = as.vector(q)
    )
  })
  do.call(rbind, rows)
}

#' Net horizontal inflow of every box per layer
#'
#' Exchanges with box -1 (the boundary of CROCO) count as inflow or outflow of
#' the box on the other side.
#'
#' @inheritParams box_transport
#' @return A data frame with one row per record, box and layer: `time`, `box`,
#'   `layer` and `inflow` (m3/s).
#' @export
box_net_inflow <- function(ds, grid, params, faces, bottoms, time = NULL) {
  if (is.null(time)) time <- croco_records(ds)
  nbox <- reticulate::py_to_r(faces$nbox)
  nlayer <- length(bottoms)
  rows <- lapply(time, function(t) {
    q <- transport_py(ds, grid, params, faces, bottoms, t)
    inflow <- reticulate::py_to_r(py_mod("fluxes")$net_inflow(faces, q))
    data.frame(
      time = as.integer(t),
      box = rep(seq_len(nbox) - 1L, nlayer),
      layer = rep(seq_len(nlayer), each = nbox),
      inflow = as.vector(inflow)
    )
  })
  do.call(rbind, rows)
}

#' Pairs of boxes that exchange water without a shared BGM face
#'
#' @param faces Cell faces from [cell_faces()].
#' @param bgm A BGM from [read_bgm()].
#' @return A data frame with columns `from` and `to`.
#' @export
pairs_without_face <- function(faces, bgm) {
  p <- matrix(as.integer(reticulate::py_to_r(py_mod("fluxes")$pairs_without_face(faces, bgm))), ncol = 2)
  data.frame(from = p[, 1], to = p[, 2])
}
