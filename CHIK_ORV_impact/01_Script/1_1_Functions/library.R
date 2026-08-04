# library list
# Consolidated from every library() call used across 00_Data/ and 01_Script/.
# source("01_Script/1_1_Functions/library.R") at the top of a script to load everything.

# ---- data wrangling ----
library(dplyr)
library(tidyr)
library(tibble)
library(purrr)
library(stringr)
library(forcats)
library(lubridate)
library(rlang)
library(glue)
library(data.table)
library(vctrs)
library(varhandle)
library(class)

# ---- I/O ----
library(readxl)
library(writexl)
library(reactable)
library(htmlwidgets)

# ---- visualization ----
library(ggplot2)
library(scales)
library(patchwork)
library(cowplot)
library(ggpubr)
library(ggsci)
library(ggridges)
library(ggnewscale)
library(corrplot)
library(scico)
library(paletteer)
library(RColorBrewer)
library(viridis)
library(viridisLite)
library(colorspace)
library(nord)
library(hrbrthemes)
library(vioplot)
library(png)
library(tiff)
library(grid)

# ---- spatial ----
library(sf)
library(raster)
library(exactextractr)
library(fasterize)
library(geodata)
library(tidyterra)
library(tmap)
library(maps)
library(rnaturalearth)
library(rnaturalearthdata)
library(ncdf4)
library(fields)

# ---- statistics ----
library(mgcv)
library(car)
library(betareg)
library(statmod)
library(MASS)
library(truncnorm)
library(FNN)
library(NHANES)

# ---- machine learning ----
library(randomForest)
library(gbm)
library(ranger)
library(caret)
library(tidymodels)
library(parsnip)
library(recipes)
library(workflows)
library(modeldata)
library(yardstick)
library(vip)
library(h2o)
library(blockCV)
library(spatialRF)
library(spatialsample)
library(ecospat)

# ---- bayesian ----
library(rstan)
library(rstanarm)
library(bayesplot)

# ---- ODE / simulation ----
library(deSolve)
library(lhs)

# ---- parallel computing ----
library(parallel)
library(doParallel)
library(snowfall)

# ---- misc ----
library(reshape)

options(scipen = 999)
