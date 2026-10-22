# croco2atlantis para R

Interfaz de R al paquete de Python `croco2atlantis`. Todo el cálculo se hace en Python mediante [reticulate](https://rstudio.github.io/reticulate/); el paquete de R solo lo llama y devuelve los resultados como data frames, matrices u objetos `sf`. El conjunto de datos de CROCO, la malla, el BGM y la asignación de celdas se quedan como objetos de Python, así que no se copian arreglos grandes entre los dos lenguajes.

## Instalación

```r
# paquete de R, desde GitHub
remotes::install_github("rcavieses/CROCOtoATLANTIS", subdir = "r")

# paquete de Python, en un virtualenv llamado "croco2atlantis" (una sola vez)
croco2atlantis::install_croco2atlantis()
```

Desde una copia local del repositorio: `R CMD INSTALL r` y `croco2atlantis::install_croco2atlantis(source = ".", editable = TRUE)`.

Al cargarse, el paquete usa el virtualenv `croco2atlantis` si existe (otro nombre con la variable de entorno `CROCO2ATLANTIS_VENV`), salvo que `RETICULATE_PYTHON` esté definida. Otro Python: `use_croco2atlantis(python = "/ruta/a/python")` antes de cualquier otra llamada.

## Uso

```r
library(croco2atlantis)

bottoms <- c(18, 54, 270, 684, 918, 1782, 3100)
bgm <- read_bgm("modelo.bgm")
ds <- croco_open("croco_avg.nc")
grid <- croco_grid(ds)
params <- croco_vertical_params(ds)

boxes <- bgm_boxes(bgm, bottoms)            # sf: cajas con profundidad y número de capas
cellmap <- assign_cells(grid, bgm)          # una vez
faces <- cell_faces(cellmap, grid)          # una vez

means <- box_layer_means(ds, grid, params, cellmap, bgm, bottoms,
                         variables = c("temp", "salt"), time = 1:2)
q <- box_transport(ds, grid, params, faces, bottoms, time = 1)
inflow <- box_net_inflow(ds, grid, params, faces, bottoms, time = 1)
fechas <- croco_times(ds, origin = "2000-01-01")
```

## Convenciones

- `box` es el número de caja de ATLANTIS, desde 0 como en el BGM. La caja −1 es el exterior de las cajas (las fronteras abiertas de CROCO).
- `layer` es la posición de la capa en `bottoms`, desde 1 en la superficie. En los archivos de ATLANTIS la numeración va al revés (la capa 0 es la del fondo); ese cambio lo hará el escritor.
- `time` es el registro del archivo de CROCO, desde 1.
- `box_layer_means` devuelve solo las capas que tiene cada caja; `filled` indica los valores copiados de otra caja y `volume` el volumen de agua de CROCO promediado.

## Funciones

| Función | Qué devuelve |
|---|---|
| `install_croco2atlantis`, `use_croco2atlantis`, `croco2atlantis_available`, `croco2atlantis_version` | Instalación y elección del Python |
| `read_bgm`, `bgm_boxes`, `bgm_faces` | BGM (objeto de Python) y sus cajas y caras como `sf` |
| `croco_open`, `croco_grid`, `croco_vertical_params`, `croco_times`, `croco_s_depths` | Archivo de CROCO, malla, parámetros verticales, fechas y profundidad de los niveles s |
| `assign_cells`, `cell_boxes` | Asignación de celdas a cajas y su matriz en R |
| `box_layer_means` | Medias por caja y capa ponderadas por volumen, con relleno opcional |
| `cell_faces`, `box_transport`, `box_net_inflow`, `pairs_without_face` | Transporte entre cajas, balance por caja y pares sin cara en el BGM |

La ayuda de cada función (`?box_layer_means`) explica los argumentos; el método está en el [README](../README.md) del proyecto.

## Pruebas

```r
testthat::test_local("r")
```

Las pruebas usan un BGM y un archivo de CROCO sintéticos y se saltan si Python no tiene `croco2atlantis`.
