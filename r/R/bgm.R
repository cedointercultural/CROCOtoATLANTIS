# ATLANTIS geometry (BGM).

#' Read an ATLANTIS BGM file
#'
#' @param path Path of the BGM file.
#' @return The Python `Bgm` object, to pass to the other functions. Use
#'   [bgm_boxes()] and [bgm_faces()] to get its boxes and faces in R.
#' @export
read_bgm <- function(path) {
  py_mod("geometry.bgm")$read_bgm(normalizePath(path, mustWork = TRUE))
}

bgm_crs <- function(bgm) sf::st_crs(reticulate::py_to_r(bgm$projection))

#' Boxes of a BGM as an sf object
#'
#' @param bgm A BGM from [read_bgm()].
#' @param bottoms Optional layer bottoms (m); when given, adds the number of
#'   layers of every box (`n_layers`).
#' @return An `sf` data frame with one row per box: `box` (id from 0), `label`,
#'   `depth` (mean depth, m, positive; 0 for islands), `area` (m2), `vertmix`,
#'   `horizmix`, `inside_x`, `inside_y` (interior point) and the polygon, in the
#'   projection of the BGM.
#' @export
bgm_boxes <- function(bgm, bottoms = NULL) {
  boxes <- reticulate::py_to_r(bgm$boxes)
  get <- function(attr) vapply(boxes, function(b) reticulate::py_to_r(reticulate::py_get_attr(b, attr)), numeric(1))
  inside <- t(vapply(boxes, function(b) as.numeric(reticulate::py_to_r(b$inside)), numeric(2)))
  polygons <- lapply(boxes, function(b) sf::st_polygon(list(reticulate::py_to_r(b$vert))))
  out <- data.frame(
    box = as.integer(get("id")),
    label = vapply(boxes, function(b) reticulate::py_to_r(b$label), character(1)),
    depth = as.numeric(reticulate::py_to_r(bgm$depth)),
    area = get("area"),
    vertmix = get("vertmix"),
    horizmix = get("horizmix"),
    inside_x = inside[, 1],
    inside_y = inside[, 2]
  )
  if (!is.null(bottoms)) {
    out$n_layers <- as.integer(reticulate::py_to_r(
      py_mod("vertical")$n_layers(bgm$depth, np_float(bottoms))
    ))
  }
  sf::st_sf(out, geometry = sf::st_sfc(polygons, crs = bgm_crs(bgm)))
}

#' Faces of a BGM as an sf object
#'
#' @param bgm A BGM from [read_bgm()].
#' @return An `sf` data frame with one row per face: `face` (id from 0), `left`
#'   and `right` (the boxes on each side looking from `p1` to `p2`), `length` (m)
#'   and the segment from `p1` to `p2`.
#' @export
bgm_faces <- function(bgm) {
  faces <- reticulate::py_to_r(bgm$faces)
  lr <- matrix(reticulate::py_to_r(bgm$face_lr()), ncol = 2)
  lines <- lapply(faces, function(f) {
    sf::st_linestring(rbind(reticulate::py_to_r(f$p1), reticulate::py_to_r(f$p2)))
  })
  out <- data.frame(
    face = seq_along(faces) - 1L,
    left = as.integer(lr[, 1]),
    right = as.integer(lr[, 2]),
    length = as.numeric(reticulate::py_to_r(bgm$face_lengths()))
  )
  sf::st_sf(out, geometry = sf::st_sfc(lines, crs = bgm_crs(bgm)))
}

# Interior point of every box, (n_box, 2) numpy array, used as box centroids.
box_inside <- function(bgm) {
  boxes <- reticulate::py_to_r(bgm$boxes)
  inside <- t(vapply(boxes, function(b) as.numeric(reticulate::py_to_r(b$inside)), numeric(2)))
  reticulate::np_array(inside, dtype = "float64")
}
