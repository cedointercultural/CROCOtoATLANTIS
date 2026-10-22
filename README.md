# croco2atlantis

Herramienta en Python para convertir salidas del modelo oceánico [CROCO](https://www.croco-ocean.org/) en los forzantes hidrográficos del modelo de ecosistemas [ATLANTIS](https://github.com/Atlantis-Ecosystem-Model): temperatura, salinidad e intercambio de agua entre cajas y capas.

**Autor:** Ricardo Cavieses.

## Estado

**En desarrollo.** La conversión completa todavía no funciona y la CLI solo muestra la versión.

| Módulo | Estado | Qué hace |
|---|---|---|
| `readers/croco.py` | Hecho y probado | Abre un archivo de CROCO, lee la malla (con el anillo fantasma enmascarado y la frontera abierta aparte), los parámetros verticales y las fechas, y calcula la profundidad de los niveles s (`Vtransform` 1 y 2) |
| `geometry/bgm.py` | Hecho y probado | Lee el BGM: cajas, caras, normales, proyección |
| `vertical.py` | Hecho y probado | Capas de cada caja, solapamiento de niveles s con capas y remapeo ponderado por espesor |
| `horizontal.py` | Hecho y probado | Asigna celdas a cajas, calcula medias por caja y capa ponderadas por volumen y rellena las cajas sin datos |
| `fluxes.py` | Hecho y probado | Transporte de agua entre pares de cajas y por capa, en las caras de las celdas de CROCO |
| `exchanges.py` | Pendiente | Origen y destino de cada intercambio, flujo vertical, corrección de hiperdifusión, paso a m³ por paso de tiempo y orden de capas de ATLANTIS |
| `writers/atlantis_nc.py` | Pendiente | Escritura de `exchange`, `dest_b`, `dest_k`, temperatura y salinidad |
| `validate.py` | Pendiente | Volumen, continuidad, rangos físicos y orden de capas |
| `cli.py` | Pendiente | CLI con configuración en YAML |

## Instalación y pruebas

```
pip install -e ".[dev]"
pytest
```

Sin instalar: `PYTHONPATH=src python3 -m pytest`. Una sola prueba: `pytest tests/test_fluxes.py::test_nondivergent_flow_closes_every_box`.

## Problema

CROCO usa niveles verticales que siguen el fondo (coordenadas s), con una malla curvilínea. ATLANTIS usa polígonos (cajas) con capas de espesor fijo y caras compartidas entre cajas. Los valores no se pueden copiar de una malla a otra: hay que reasignarlos conservando el contenido de calor y sal y el volumen de agua.

## Método

### Datos de entrada: salidas promediadas (`avg`)

La herramienta usa los archivos de promedios de CROCO (`avg`), no los instantáneos (`his`). En una corrida con marea, una foto cada 12 h submuestrea la marea: M2 (12.42 h) aparece como una oscilación falsa de unos 15 días, y S2 (12 h exactas) se muestrea siempre en la misma fase y queda como un sesgo constante del nivel del mar. Los dos tipos de archivo tienen la misma estructura, así que el lector sirve para ambos. En los `avg`, el tiempo de cada registro es el centro de su ventana de promedio.

Un promedio de 12 h elimina casi del todo las mareas semidiurnas, pero deja buena parte de las diurnas (K1, O1).

### Vertical

- Las capas de ATLANTIS se definen por la profundidad de su fondo, por ejemplo `(18, 54, 270, 684, 918, 1782, 3100)` m.
- Cada caja conserva las capas cuyo tope está por encima de su profundidad media (`botz` del BGM), y la última se trunca en esa profundidad.
- Para promediar, la capa más profunda de cada caja toma toda la columna de agua bajo su tope (opción `deepest_to_bottom`), porque algunas columnas del modelo son más profundas que la media de la caja.
- La primera capa empieza en la superficie libre, así que incluye la elevación `zeta`.
- Los niveles de CROCO se reparten entre las capas según su solapamiento real en profundidad.
- Convención interna: en CROCO el índice 0 es el fondo; en las capas de ATLANTIS el índice 0 es la más somera hasta que se escriben los archivos, donde la capa 0 pasa a ser el fondo de cada caja.

### Horizontal

- Cada celda de agua de CROCO va a la caja que contiene su centro. Si el centro cae justo en el borde entre dos cajas, gana la primera.
- Las celdas de agua cuyo centro queda fuera de todas las cajas, o dentro de una isla del BGM que la malla de CROCO no resuelve, pasan a la caja con agua más cercana. Así no se pierde agua ni flujo. Un `max_distance` opcional deja sin caja las que estén más lejos.
- El valor de una caja en una capa es la media de sus celdas ponderada por volumen (área de la celda × espesor de su agua dentro de la capa).
- Un valor ausente (NaN) en una celda de agua es un error y detiene el cálculo. En CROCO la tierra vale 0, no NaN, por eso siempre se enmascara con `mask_rho`.
- Las cajas o capas sin datos (cajas sin ninguna celda, o capas más profundas que todas sus columnas) se rellenan con el perfil de la caja más cercana que tenga datos en **todas** las capas que faltan. Así la columna rellenada es un perfil real y no capas de cajas distintas. Un `max_distance` opcional evita copiar valores de muy lejos, por ejemplo en cajas fuera del dominio de CROCO.

### Transporte entre cajas

- El transporte se calcula en las caras de las celdas de la malla C de CROCO, donde viven `u` y `v`. Cada cara entre dos celdas de cajas distintas aporta velocidad × longitud de la cara × espesor del agua, integrado en cada capa. Los niveles de la cara son la media de los de sus dos celdas, como en CROCO.
- Así, el balance de agua de cada caja es la divergencia del propio modelo y es coherente con los volúmenes usados en los promedios.
- No se muestrea la velocidad a lo largo de las caras del BGM: probado con datos reales, ese método produce en el balance neto de cada caja un error del mismo orden que el balance real.
- El resultado es por par de cajas (`pairs`), con signo positivo de `pairs[:, 0]` a `pairs[:, 1]`.
- La caja −1 representa el exterior de las cajas: las fronteras abiertas de CROCO (celdas de agua del anillo fantasma) y cualquier celda de agua sin caja.
- Las capas del transporte son las estándar, de la superficie al fondo, sin truncar. El paso a las capas de cada caja se hace al construir los intercambios.
- Con archivos `avg`, el transporte es la velocidad media por el espesor del `zeta` medio. Omite la correlación entre corriente y elevación de marea (transporte de Stokes de marea), que CROCO no guarda.

### Uso actual (API)

```python
from croco2atlantis.readers import croco
from croco2atlantis.geometry.bgm import read_bgm
from croco2atlantis import horizontal as hz, fluxes

bottoms = [18, 54, 270, 684, 918, 1782, 3100]
ds = croco.open_output("croco_avg.nc")
grid = croco.read_grid(ds)
params = croco.read_vertical_params(ds)
bgm = read_bgm("modelo.bgm")

cellmap = hz.assign_cells_from_grid(grid, bgm)   # una vez
faces = fluxes.cell_faces(cellmap, grid)         # una vez

t = 0
z_w = croco.s_depths(grid.h, ds.zeta.isel(time=t).values, params, "w")
weights = hz.box_layer_weights(z_w, cellmap, grid.cell_area, bgm.depth, bottoms)
temp = weights.apply(ds.temp.isel(time=t).values)            # (caja, capa)
temp, filled = hz.fill_missing(temp, bgm.depth, bottoms, centroids)
q = fluxes.face_transport(faces, ds.u.isel(time=t).values, ds.v.isel(time=t).values, z_w, bottoms)
inflow = fluxes.net_inflow(faces, q)                          # balance por caja y capa
```

`centroids` es un arreglo (caja, 2) en metros en la proyección del BGM, por ejemplo los puntos interiores de las cajas: `numpy.array([b.inside for b in bgm.boxes])`.

## Decisiones pendientes

- **Frontera abierta:** a qué caja de ATLANTIS se asigna el agua que entra o sale por las fronteras abiertas de CROCO (caja −1): la caja vecina fuera del dominio o una caja de frontera.
- **Pares sin cara en el BGM:** como las celdas se asignan por su centro, dos cajas que solo se tocan en un vértice, o que separa una caja delgada, pueden intercambiar agua sin compartir una cara (`fluxes.pairs_without_face`). Falta decidir si quedan como intercambio directo o se reparten por la caja intermedia. La corrección de hiperdifusión necesita la geometría de una cara.
- **Mareas diurnas:** si se promedia además a 24 h o más (por ejemplo, la media móvil de dos registros de 12 h) antes de construir los flujos.
- **Flujo vertical:** por continuidad a partir de los flujos horizontales o a partir de `w`.
- **Corrección de hiperdifusión:** dividir por la distancia a través de la caja según el ángulo de la cara, por la longitud o por el área.
- **Cajas fuera del dominio de CROCO:** si se rellenan (y con qué `max_distance`) o se dejan sin datos.

## Registro de cambios

### 2026-10-22: revisión del código y la lógica

- Datos de entrada: se pasa de `his` a `avg` por el submuestreo de la marea. `open_history` se renombró a `open_output`.
- Transporte: se reemplazó el muestreo de la velocidad de la celda más cercana a lo largo de las caras del BGM por el transporte en las caras de las celdas de la malla C, que cierra el balance de cada caja. Se agregaron `cell_faces`, `net_inflow` y `pairs_without_face`, y se eliminaron `rho_velocity` y `face_samples`.
- Frontera abierta: `CrocoGrid.open_boundary` guarda las celdas de agua del anillo fantasma; su transporte va a la caja −1.
- Celdas sin caja o en islas del BGM: ahora pasan a la caja con agua más cercana (`move_to_water_boxes`); antes se descartaban.
- Relleno: cada caja toma todas sus capas faltantes de una sola caja; antes cada capa podía venir de una caja distinta.
- Asignación de celdas: un centro sobre el borde entre dos cajas va a la primera; antes no iba a ninguna.
- `BoxLayerWeights.apply` acepta cualquier arreglo y falla si hay NaN en el agua; `s_depths` solo acepta `"w"` o `"rho"`.

## Alcance

- Solo CROCO; no pretende cubrir otros modelos ROMS.
- Pensado para que otros investigadores que usen CROCO lo reutilicen con su propio dominio y su propia geometría de ATLANTIS (archivo BGM). Todo lo que depende de una región (BGM, espesores de capa, nombres de archivo) irá en un archivo de configuración YAML, no en el código.
- No usa Hydro ni HydroConstruct. El método sigue el de proyectos existentes desde ROMS y otros modelos (NEMOW Climate Data Workshop, ROMS_Atlantis_translation, SSMtoAtlantis).
- Uso interno por ahora; el código se publica en GitHub.

## Estructura

```
src/croco2atlantis/
  readers/croco.py      # único módulo que conoce CROCO
  geometry/bgm.py       # lectura del BGM: cajas, caras, normales
  vertical.py           # capas de ATLANTIS y remapeo vertical
  horizontal.py         # celda→caja, promedio por caja y relleno
  fluxes.py             # transporte entre cajas en las caras de las celdas
  exchanges.py          # origen/destino, flujo vertical, hiperdifusión y orden de capas
  writers/atlantis_nc.py
  validate.py
  cli.py
tests/
examples/
```

## Licencia

GPL-3.0 (ver [LICENSE](LICENSE)).
