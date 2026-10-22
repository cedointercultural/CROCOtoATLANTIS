# Access to the Python package.

.c2a <- new.env(parent = emptyenv())

default_envname <- function() Sys.getenv("CROCO2ATLANTIS_VENV", "croco2atlantis")

.onLoad <- function(libname, pkgname) {
  if (Sys.getenv("RETICULATE_PYTHON") == "" && reticulate::virtualenv_exists(default_envname())) {
    reticulate::use_virtualenv(default_envname(), required = FALSE)
  }
}

#' Install the Python package in a virtualenv
#'
#' Creates the virtualenv `envname` if it does not exist and installs
#' `croco2atlantis` with pip. Restart R afterwards if Python was already
#' running with another environment.
#'
#' @param envname Name or path of the virtualenv.
#' @param source What pip installs: the GitHub repository by default, or the
#'   path of a local copy of the repository.
#' @param editable Install `source` in editable mode (`pip install -e`), for a
#'   local copy that is being developed.
#' @param python Python used to create the virtualenv, `NULL` for the default.
#' @return The path of the Python of the virtualenv, invisibly.
#' @export
install_croco2atlantis <- function(envname = default_envname(),
                                   source = "git+https://github.com/rcavieses/CROCOtoATLANTIS",
                                   editable = FALSE, python = NULL) {
  if (!reticulate::virtualenv_exists(envname)) {
    reticulate::virtualenv_create(envname, python = python)
  }
  py <- reticulate::virtualenv_python(envname)
  if (!grepl("^git\\+", source) && file.exists(source)) {
    source <- normalizePath(source)
  }
  args <- c("-m", "pip", "install", if (editable) "-e", shQuote(source))
  status <- system2(py, args)
  if (status != 0) stop("pip install failed with status ", status)
  invisible(py)
}

#' Choose the Python that has croco2atlantis
#'
#' Must be called before the first call to Python in the session.
#'
#' @param envname Name or path of a virtualenv, used when `python` is `NULL`.
#' @param python Path of a Python executable.
#' @return The Python configuration from [reticulate::py_config()], invisibly.
#' @export
use_croco2atlantis <- function(envname = default_envname(), python = NULL) {
  if (is.null(python)) {
    reticulate::use_virtualenv(envname, required = TRUE)
  } else {
    reticulate::use_python(python, required = TRUE)
  }
  invisible(reticulate::py_config())
}

#' Is the Python package available?
#'
#' @return `TRUE` when the Python in use can import `croco2atlantis`.
#' @export
croco2atlantis_available <- function() {
  reticulate::py_module_available("croco2atlantis")
}

#' Version of the Python package
#'
#' @return A string, for example `"0.0.1"`.
#' @export
croco2atlantis_version <- function() {
  reticulate::py_to_r(py_mod()$`__version__`)
}

# Python module croco2atlantis or one of its submodules, imported once without
# conversion so that the objects stay in Python.
py_mod <- function(name = NULL) {
  full <- if (is.null(name)) "croco2atlantis" else paste0("croco2atlantis.", name)
  if (is.null(.c2a[[full]])) {
    if (!croco2atlantis_available()) {
      stop("the Python package croco2atlantis is not available in ",
           reticulate::py_exe(), "; install it with install_croco2atlantis()",
           " or choose its Python with use_croco2atlantis()", call. = FALSE)
    }
    .c2a[[full]] <- reticulate::import(full, convert = FALSE)
  }
  .c2a[[full]]
}

np_float <- function(x) reticulate::np_array(as.numeric(x), dtype = "float64")

# R record index (from 1) to Python index (from 0).
py_index <- function(time) {
  if (length(time) != 1 || is.na(time) || time < 1 || time != round(time)) {
    stop("time must be a record index from 1, found ", format(time), call. = FALSE)
  }
  as.integer(time) - 1L
}

null_or_float <- function(x) if (is.null(x)) NULL else as.numeric(x)
