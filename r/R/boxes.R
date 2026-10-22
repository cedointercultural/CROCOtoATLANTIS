# CROCO cells to ATLANTIS boxes, and box-layer means.

#' Assign the CROCO cells to the boxes of a BGM
#'
#' Every water cell goes to the box that contains its centre; water cells
#' outside every box or inside an island of the BGM go to the nearest box with
#' water. Call it once per grid and BGM.
#'
#' @param grid A grid from [croco_grid()].
#' @param bgm A BGM from [read_bgm()].
#' @param max_distance Optional distance (m): cells farther than that from every
#'   box with water are left without a box.
#' @return The Python `CellBoxMap` object. Use [cell_boxes()] to see it in R.
#' @export
assign_cells <- function(grid, bgm, max_distance = NULL) {
  py_mod("horizontal")$assign_cells_from_grid(grid, bgm, max_distance = null_or_float(max_distance))
}

#' Box of every CROCO cell
#'
#' @param cellmap A cell map from [assign_cells()].
#' @return An integer matrix (eta, xi) with the box of every cell, -1 for land,
#'   the ghost ring and cells without a box. Attribute `moved` is a logical
#'   matrix marking the water cells moved to the nearest box with water.
#' @export
cell_boxes <- function(cellmap) {
  box <- reticulate::py_to_r(cellmap$box)
  storage.mode(box) <- "integer"
  moved <- reticulate::py_to_r(cellmap$moved)
  if (!is.null(moved)) attr(box, "moved") <- moved
  box
}

#' Volume-weighted means of CROCO variables per box and layer
#'
#' For every record, the mean of each variable over the cells of each box,
#' weighted by the volume of every cell inside each layer. Box-layers without
#' data can be filled with the profile of the nearest box that has them.
#'
#' @param ds A dataset from [croco_open()].
#' @param grid A grid from [croco_grid()].
#' @param params Parameters from [croco_vertical_params()].
#' @param cellmap A cell map from [assign_cells()].
#' @param bgm A BGM from [read_bgm()].
#' @param bottoms Depth of the bottom of every layer, m, positive.
#' @param variables Names of the CROCO variables on rho points (time, s_rho,
#'   eta_rho, xi_rho), such as `"temp"` and `"salt"`.
#' @param time Records, from 1; `NULL` for all.
#' @param fill Fill the box-layers without data (see Details).
#' @param max_distance With `fill`, do not fill from a box farther than this
#'   (m), so that boxes outside the CROCO domain stay `NA`.
#' @param deepest_to_bottom Average the whole water column below the top of the
#'   deepest layer of every box, not just down to the mean depth of the box.
#' @details All the missing layers of a box are filled from a single box: the
#'   nearest one (by the interior points of the BGM) with data in all of them.
#' @return A data frame with one row per record, box, layer and variable, only
#'   for the layers each box has: `time`, `box`, `layer`, `variable`, `value`
#'   (`NA` where there is no data), `filled` (the value was copied from another
#'   box) and `volume` (m3 of CROCO water averaged, 0 for filled values).
#' @export
box_layer_means <- function(ds, grid, params, cellmap, bgm, bottoms,
                            variables = c("temp", "salt"), time = NULL, fill = TRUE,
                            max_distance = NULL, deepest_to_bottom = TRUE) {
  hz <- py_mod("horizontal")
  if (is.null(time)) time <- croco_records(ds)
  b <- np_float(bottoms)
  depth <- bgm$depth
  nbox <- reticulate::py_to_r(bgm$nbox)
  nlayer <- length(bottoms)
  n_layers <- reticulate::py_to_r(py_mod("vertical")$n_layers(depth, b))
  index <- expand.grid(box = seq_len(nbox) - 1L, layer = seq_len(nlayer))
  keep <- index$layer <= n_layers[index$box + 1L]
  centroids <- if (fill) box_inside(bgm)

  rows <- list()
  for (t in time) {
    weights <- hz$box_layer_weights(
      s_depths_py(ds, grid, params, t), cellmap, grid$cell_area, depth, b,
      deepest_to_bottom = deepest_to_bottom
    )
    volume <- reticulate::py_to_r(weights$volume)
    for (v in variables) {
      means <- weights$apply(croco_field(ds, v, t))
      filled <- matrix(FALSE, nbox, nlayer)
      if (fill) {
        res <- hz$fill_missing(means, depth, b, centroids, max_distance = null_or_float(max_distance))
        means <- reticulate::py_get_item(res, 0L)
        filled <- reticulate::py_to_r(reticulate::py_get_item(res, 1L))
      }
      rows[[length(rows) + 1L]] <- data.frame(
        time = as.integer(t),
        index[keep, ],
        variable = v,
        value = as.vector(reticulate::py_to_r(means))[keep],
        filled = as.vector(filled)[keep],
        volume = as.vector(volume)[keep]
      )
    }
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}
