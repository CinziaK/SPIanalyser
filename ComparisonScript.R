AnalyseSPIData = function(GBP.file,GOI.file, k = 0.4, writename = "Merge.file.csv"){
  # GOIfile and GBP file should be in form ".txt"


GBP.file = read_csv(GBP.file)
GOI.file = read_csv(GOI.file) 


Merge.file <- inner_join(GBP.file, GOI.file, by = c("Plate", "Row", "Column", "ORF"), suffix = c(".GBP",".GOI"))
Merge.file 


Mean_LGR = c() # calculates the mean LGR of both comparisons 
Mean_LGR = rowMeans(Merge.file[c('mean_LGR.GBP', 'mean_LGR.GOI')], na.rm=TRUE)
Mean_Z = c() # calculates the mean Z score of both comparisons
Mean_Z = rowMeans(Merge.file[c('mean_Z_score.GBP', 'mean_Z_score.GOI')], na.rm=TRUE)

Merge.File = mutate( Merge.file, Mean_LGR, Mean_Z)
Merge.File #creative

write.csv(Merge.File, file = sprintf("Merge.File_%s.csv",Var))


SPIplots = function(Merge.File){
  A = arrange(Merge.File, desc(Mean_LGR)) #arrange largest to smallest
  lim = ceiling(max(A$Mean_LGR)) # find best limits for graph
  plot(A$Mean_LGR, xlab = "GFP strains", ylab = "LGR", pch = 19,col = rainbow(1),cex = .4,panel.first=grid(20,10),ylim = c(-1*lim,lim))
  abline(h = c(k,-k)) # add lines at ±2
  
  
  B = filter(Merge.File, (mean_LGR.GOI != "n/a")&(mean_LGR.GBP != "n/a"))
  Xlim = ceiling(as.numeric(max(B$mean_LGR.GOI))) # Best limits for graph
  Ylim = ceiling(as.numeric(max(B$mean_LGR.GBP)))
  Color = (9*(B$Mean_LGR>= k)+1) # Vector of colors, red = SPI, black  = non-SPI
  Color = 9*(B$Mean_LGR<=- k)+1
  Color = (9*(B$Mean_LGR>= k)+1) 
  plot(B$mean_LGR.GOI ,B$mean_LGR.GBP ,pch = 19,panel.first=grid(12,12),asp =1,cex = .4, xlim = c(-1*Xlim,Xlim),ylim = c(-1*Ylim,Ylim), xlab = "GOI control LGR", ylab = "GBP control LGR", col = Color, main = "Control Correlation" )
  abline(h = k ,v= k)
  lines(c((k*102),(-k*100)),c((-k*100),(k*102))) # Line showing (x+y)/2 = 2, line showing SPI limit
  legend("bottomleft",c("SPI","non-SPI"), col = c(10,1), pch =19)
  
  
}
SPIplots(Merge.File) 



return(Merge.File)

}
