Readme: This file helps you to use the SPIdataShellScript.R package.

You should create a folder containing the files included in the zipped archive, these are:

2ndRCSmoother-16x.plx
areas2LGRs.R
ColAreas_extractor-4.pl
ComparisonScript.R
keyfile_letters_to_numbers.pl
SPIdataShellScipt.R
keyfile.txt

You should replace the "keyfile.txt" with your own keyfile. A keyfile is a table listing the postion of ORF identifiers associated with each position on each plate.
The keyfile should look like this...

Plate #	Row	Column	ORF
1		A	1		YML032C
etc.....

or like this
Plate #	Row	Column	ORF
1		1	1		YML032C
etc.

You should copy your "colonyAreas.txt" file into the same folder (this file should be the output from cm_engine).

To run the R script open an program that can run R scripts! We use "R Studio" on mac computers. If packages can't be installed in R Studio, try install the package in R (outside of R Studio) or visit https://support.rstudio.com/hc/en-us/articles/200554786-Problem-Installing-Packages. 

Drag the "SPIdataShellScipt.R" file from the folder onto the R Studio window, the following text should appear in the active window

>>
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


######### combining GOI and GBP control in a merge.file

source("ComparisonScript.R")

GBP.file = "mean_file_GBP.csv" # enter name of file that you want to compare 
GOI.file = "mean_file_GOI.csv" # enter name of file that you want to compare 
Var = "Fusion" # Variable added to output file name, resulting in Merge.File_"...".csv
H = AnalyseSPIData(GBP.file, GOI.file, k = 0.4)
# select cutoff value for SPIs in graph (k)
# results in SPI plot and control comparison plot 

>>

Each time you restart R studio you should run each line of code starting at the top with "install.packages...etc"
To do this in R studio, you position the cursor on the line of text, and press Control + Enter
