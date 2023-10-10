Análisis de datos de Micobacterias ambientales presentes en la Planta de
tratamiento de aguas residuales Palandacocha
================
Silvana Morocho, Jefferson Pastuña
2023-10-09

- <a href="#introducción" id="toc-introducción">Introducción</a>
- <a href="#identificación-molecular"
  id="toc-identificación-molecular">Identificación Molecular</a>
- <a href="#estadística" id="toc-estadística">Estadística</a>
  - <a href="#diversidad-alfa" id="toc-diversidad-alfa">Diversidad Alfa</a>
    - <a href="#índice-de-simpson" id="toc-índice-de-simpson">Índice de
      Simpson</a>
    - <a href="#índice-de-shannon" id="toc-índice-de-shannon">Índice de
      Shannon</a>
  - <a href="#diversidad-beta" id="toc-diversidad-beta">Diversidad Beta</a>

# Introducción

General…

# Identificación Molecular

You can include R code in the document as follows:

``` r
summary(cars)
```

    ##      speed           dist       
    ##  Min.   : 4.0   Min.   :  2.00  
    ##  1st Qu.:12.0   1st Qu.: 26.00  
    ##  Median :15.0   Median : 36.00  
    ##  Mean   :15.4   Mean   : 42.98  
    ##  3rd Qu.:19.0   3rd Qu.: 56.00  
    ##  Max.   :25.0   Max.   :120.00

# Estadística

Se usarán índices de diversidad en cada punto de muestreo para un
seguimiento al tratamiento de las aguas residuales que realiza la planta
Palandacocha. En cambio, todas las especies identificadas en la planta
serán consideradas como la diversidad de Micobacterias presentes en las
aguas residuales del cantón Tena (representado por aquellas que llegan a
la planta Palandacocha).

## Diversidad Alfa

Los índices de diversidad más usados en la literatura son Simpson y
Shannon (cita). El índice de Simpson está relacionado con la dominancia
de las especies en un área o ecosistema en particular (cita). En este
caso se calculará el índice de Simpson en muestras (área o punto de
muestreo) antes del ingreso, dentro y a la salida de la planta de
tratamiento Palandacocha. De existir una especie dominante a la salida
de la planta podría deberse a la resistencia del microorganismo a los
métodos descontaminantes empleadas en la planta.

Por otro lado, el índice de Shannon está relacionado a la
equidad/uniformidad de las especies en un área o ecosistema en
particular (cita). Este estadístico nos permitirá conocer si las
diferentes áreas muestrales son iguales, más o menos uniformes. De esta
manera, la uniformidad semejante entre los sitios muéstreles podría
deberse a la resistencia de los microrganismos a los métodos de
descontaminación empleadas en la planta de Palandacocha.

Para el cálculo de la diversidad alfa se usará el paquete R
vegan(<https://github.com/vegandevs/vegan/tree/master>) disponible en
GitHub (vegandevs/vegan). En la siguiente línea de código se instalará y
activará las bibliotecas del paquete R vegan.

### Índice de Simpson

Descripción

### Índice de Shannon

Descripción

## Diversidad Beta

Descripción

![](PTAR_Palandacocha_files/figure-gfm/pressure-1.png)<!-- -->

Note that the `echo = FALSE` parameter was added to the code chunk to
prevent printing of the R code that generated the plot.
