###### combining GOI and GBP control in a merge.file. For comparison of 2 output files from areas2LGR.R script.

source("ComparisonScript.R")

GBP.file = "mean_file_GBP.csv" # enter name of file that you want to compare 
GOI.file = "mean_file_GOI.csv" # enter name of file that you want to compare 
Var = "Fusion" # Variable added to output file name, resulting in Merge.File_"...".csv
H = AnalyseSPIData(GBP.file, GOI.file, k = 0.4)
# select cutoff value for SPIs in graph (k)
# results in SPI plot and control comparison plot 
