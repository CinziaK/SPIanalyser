#Shell script for SPI analysis. Works together with areas2LGRs.R, ColAreas_extractor-4.pl & 2ndRCsmoother-16x.pl

install.packages("tidyverse") # make sure to install package before first use
library(tidyverse) # load library


setwd("/Users/thorpelab/Desktop/rscript") # set working directory
source("areas2LGRs.R") # load function

query = "fusiongbp" # name of query strain, needs to be identical as in colonyAreas.txt (e.g.pHT...)
control = "pHT4" # name of control strain, needs to be identical as in colonyAreas.txt (e.g.pHT...)
Var = "GBP" # Variable added to output file name, results in mean_file_"...".csv and replicates_file_"...".csv

areas2LGRs("colonyAreas.txt", query, control, smooth = T, write = TRUE)
# optional smoothing(T/F): for screens where the majority of colonies are not inhibited in growth, corrects for non-specific plate effects, unsmoothed data still in output file 


