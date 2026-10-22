setup <- function() {
  dir <- withr_tempdir()
  bgm <- read_bgm(synthetic_bgm(dir))
  ds <- croco_open(synthetic_croco(dir))
  grid <- croco_grid(ds)
  list(bgm = bgm, ds = ds, grid = grid, params = croco_vertical_params(ds),
       cellmap = assign_cells(grid, bgm))
}

withr_tempdir <- function() {
  dir <- tempfile("c2a")
  dir.create(dir)
  dir
}

test_that("the Python package is found", {
  skip_if_no_python()
  expect_match(croco2atlantis_version(), "^[0-9]+\\.[0-9]+")
})

test_that("BGM boxes and faces come back as sf", {
  skip_if_no_python()
  bgm <- read_bgm(synthetic_bgm(withr_tempdir()))
  boxes <- bgm_boxes(bgm, BOTTOMS)
  expect_s3_class(boxes, "sf")
  expect_equal(boxes$box, 0:1)
  expect_equal(boxes$depth, c(100, 100))
  expect_equal(boxes$n_layers, c(3L, 3L))
  expect_equal(as.numeric(sf::st_area(boxes)), c(2e6, 2e6))
  expect_equal(sf::st_crs(boxes), sf::st_crs(PROJ))
  faces <- bgm_faces(bgm)
  expect_equal(c(faces$left, faces$right, faces$length), c(0, 1, 2000))
})

test_that("records, depths and cells", {
  skip_if_no_python()
  s <- setup()
  expect_equal(format(croco_times(s$ds, "2000-01-01")), c("2000-01-01 00:00:00", "2000-01-01 12:00:00"))
  z <- croco_s_depths(s$ds, s$grid, s$params, time = 2)
  expect_equal(dim(z), c(5, 6, 6))
  expect_equal(z[, 3, 3], c(-100, -75, -50, -25, 0))
  expect_error(croco_s_depths(s$ds, s$grid, s$params, time = 0), "from 1")
  box <- cell_boxes(s$cellmap)
  expect_equal(dim(box), c(6, 6))
  expect_true(all(box[1, ] == -1))
  expect_true(all(box[2:5, 2:3] == 0) && all(box[2:5, 4:5] == 1))
})

test_that("box-layer means are volume weighted per box", {
  skip_if_no_python()
  s <- setup()
  m <- box_layer_means(s$ds, s$grid, s$params, s$cellmap, s$bgm, BOTTOMS)
  expect_equal(nrow(m), 2 * 2 * 3 * 2)  # records x boxes x layers x variables
  temp <- m[m$variable == "temp" & m$time == 1, ]
  expect_equal(temp$value[temp$box == 0], rep(15, 3))
  expect_equal(temp$value[temp$box == 1], rep(20, 3))
  # 8 cells of 500 m x 500 m per box; layers of 18, 36 and the remaining 46 m
  expect_equal(temp$volume[temp$box == 0], 8 * 250000 * c(18, 36, 46))
  expect_false(any(m$filled))
})

test_that("transport and net inflow between boxes", {
  skip_if_no_python()
  s <- setup()
  faces <- cell_faces(s$cellmap, s$grid)
  q <- box_transport(s$ds, s$grid, s$params, faces, BOTTOMS, time = 1)
  expect_equal(q$from, rep(0L, 3))
  expect_equal(q$to, rep(1L, 3))
  # 4 cell faces of 500 m at 1 m/s
  expect_equal(q$transport, 2000 * c(18, 36, 46))
  inflow <- box_net_inflow(s$ds, s$grid, s$params, faces, BOTTOMS, time = 1)
  expect_equal(inflow$inflow[inflow$box == 1], q$transport)
  expect_equal(inflow$inflow[inflow$box == 0], -q$transport)
  expect_equal(nrow(pairs_without_face(faces, s$bgm)), 0)
})
