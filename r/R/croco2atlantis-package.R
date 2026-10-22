#' croco2atlantis: CROCO output to ATLANTIS forcing, from R
#'
#' R interface to the `croco2atlantis` Python package. Every computation runs in
#' Python through 'reticulate'; the CROCO dataset, the grid, the BGM and the
#' cell maps stay as Python objects, and only the results come back to R as
#' data frames, matrices or 'sf' objects.
#'
#' Set up Python once with [install_croco2atlantis()]. When a virtualenv named
#' `croco2atlantis` (or the name in the `CROCO2ATLANTIS_VENV` environment
#' variable) exists, it is used when the package loads, unless `RETICULATE_PYTHON`
#' is set; otherwise choose the Python with [use_croco2atlantis()].
#'
#' @section Conventions:
#' * `box` is the ATLANTIS box id, numbered from 0 as in the BGM. Box -1 is
#'   outside the boxes: the open boundaries of CROCO and any water cell left
#'   without a box.
#' * `layer` is the position of the layer in `bottoms`, from 1 at the surface
#'   downwards. The ATLANTIS files number layers the other way (0 is the bottom
#'   of each box); that renumbering is done by the writer, not here.
#' * `time` is the record of the CROCO file, from 1.
#' * `bottoms` are the depths of the bottom of each layer, in metres, positive,
#'   for example `c(18, 54, 270, 684, 918, 1782, 3100)`.
#'
#' @keywords internal
"_PACKAGE"
