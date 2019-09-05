#This script is used for SPI data analysis and requires a colonyAreas.txt file (from cm.engine) as input file. 
#This Version is used to analyse screens with 4 replictaes per strain. For 16 replictates, use second version.
areas2LGRs = function(Areas, query, control,smooth = T, write = TRUE){
  # take in colonyAreas file from CM engine (Areas) and compare two plasmids: experiment and control, note that "keyfile.txt" must be in the same folder 
  if(Areas != "colonyAreas.txt"){ # colAreas_extractor-4.plx require the file to be called this
    data = read_delim(Areas, delim = "\t")
    write_delim(data, "colonyAreas.txt", delim = "\t")
  }
  system(paste("perl", "ColAreas_extractor-4.pl"))
  x = read_delim("results.tab", delim = "\t") #output is Colony sizes ordered by position and plasmid -> input file for R analysis
  x2 = mutate(x,`Colony size`+1) # no 0 Colonysizes, important for calculations
  x3 <- x2 [c(1:5,7)] # filters Colony size column out 
  xnospace = apply(x3[1],2,function(x)gsub('\\s+', '',x))
  x4 = x3$Plasmid <- xnospace
  x5 = transform(x3, x4)
  
  x4 = filter(x3, Plasmid %in% c(query, control)) # filters plasmids which are not used for comparison! Thus, results.tab can contain multiple plasmids 
  xrearranged <- spread (x4,Plasmid,`\`Colony size\` + 1`) #spread data to convert different plasmids from rows to columns 
  
  colnames(xrearranged) = gsub(" ","",colnames(xrearranged))
  colnames(xrearranged)[which(colnames(xrearranged)==query)] = "query"
  colnames(xrearranged)[which(colnames(xrearranged)==control)] = "control" # rename query and control plasmid!!!

  LGR_GBP = c() # necessary later (don't know why its here actually but it works)
  
  ##############################################################################
  # calculate median 
  querymedian = c()
  controlmedian = c()     #make vector functions for medians -> advantage is that table isnt affected  
  for(j in 1:12){
    tempMedian = median(as.numeric(filter(xrearranged, Plate==j, !ORF== "BLANK")$query))
    querymedian = c(querymedian, rep(tempMedian, nrow(filter(xrearranged, Plate==j))))
    tempMedian = median(as.numeric(filter(xrearranged, Plate==j, !ORF== "BLANK")$control))
    controlmedian = c(controlmedian, rep(tempMedian, nrow(filter(xrearranged, Plate==j))))
  }
  #calculate median based on different plates (filter for 1:12) and exclude blank ORFs as they would impact media, claculate for all plasmids 
  
  xmedian = mutate(xrearranged, querymedian, controlmedian)  #insert medians in new columns with mutate function
  xmedianfiltered = filter(xmedian, !ORF== "BLANK") # filter median for blank ORFs! controls are non-growing, thus median and smoothing would be incorrect if included 
  mediancorrected <- transform(xmedianfiltered, norm_query = query/querymedian,  norm_control = control/controlmedian)
  # divide coloniesizes by plate median 
  
  #################################################################################
  # claculate logarithm
  log = mutate(mediancorrected, log_norm_query = log(norm_query), log_norm_control = log(norm_control))
  
  LGR = c() 
  LGR <- (log$log_norm_control- log$log_norm_query)
  xLGR = mutate(log, LGR)
  
  ###########################################################################################
  # Smooth data with Perl script 
  if(smooth){
  xLGR_S = xLGR %>% dplyr::select(plate = Plate, row = Row, column = Column, ORF, LGR)
  write_tsv(xLGR_S, "results.txt")
  system(paste("perl", "2ndRCsmoother-16x.plx"))
  smoothed = read_tsv("smoothed.tab")
  
  smoothed = select(smoothed, Plate, Row, Column, ORF, LGR = `LGR corrected scores`) %>% mutate(Plate = as.character(Plate), Row = as.character(Row), Column = as.character(Column))
  xLGR = xLGR %>% rename(usLGR = LGR) %>% mutate(Plate = as.character(Plate), Row = as.character(Row), Column = as.character(Column))
  
  xLGR = inner_join(xLGR, smoothed, by = c("Plate", "Row", "Column", "ORF"))
  
  plot(1:nrow(xLGR), xLGR$usLGR) # plot unsmoothed and smoothed data 
  plot(1:nrow(xLGR), xLGR$LGR)
  }
  
  ################################################################
  # calculate meanRow and meanColumn to combine replicates 
  xLGR2 <- as.tibble(xLGR [,1:3])
  xLGR2 <- mutate(xLGR2, meanRow = ceiling(as.numeric(Row)/2), meanColumn = ceiling(as.numeric(Column)/2))
  xLGR3 <- inner_join(xLGR, xLGR2, by = c("Plate", "Row", "Column"))
  xLGR3 = transform(xLGR3, Plate = as.numeric(Plate))
  xLGR3 = transform(xLGR3, Plate = as.numeric(Plate))
  xLGR3 = arrange(xLGR3, Plate, meanRow, meanColumn)
  
  ###############################################################
  # organising the single replictates in rows 
  if(smooth){
  t = rep(c(1,2,3,4),times=(4331))
  }else
  {t = rep(c(1,2,3,4),times=(4339))}
  
  xt = mutate (xLGR3,t)
  
  if(smooth){
  xtreduced <- xt [c(1,4:6,9:17)]
  }else{
  xtreduced <- xt [c(1,4:6,9:16)]}
  
  if(smooth){
  xtreduced1 = gather(xtreduced,  'x', value, query, control, norm_query, norm_control, log_norm_query, log_norm_control, usLGR, LGR)
  }else{
  xtreduced1 = gather(xtreduced,  'x', value, query, control, norm_query, norm_control, log_norm_query, log_norm_control, LGR) 
  }
  
  xtreduced2 = mutate(xtreduced1, Plasmid_replicate_stats = paste(x, t))
  xtreduced3 = xtreduced2 [c(1:4,7:8)]
  y = spread (xtreduced3, Plasmid_replicate_stats, value)
  
  #################################################################
  # Z-scores for each replicate
  LGR1 = "LGR 1"
  LGR2 = "LGR 2"
  LGR3 = "LGR 3"
  LGR4 = "LGR 4"
  which(colnames(y)==LGR1)
  y1 = y[c(which(colnames(y)==LGR1),which(colnames(y)==LGR2),which(colnames(y)==LGR3),which(colnames(y)==LGR4))]
  y1mean = mean(colMeans(y1))
  y1sd = mean(apply(y1,2,sd))
  
  Z_score_1 = c()
  Z_score_2 = c()
  Z_score_3 = c()
  Z_score_4 = c()
  Z_score_1 = (y$"LGR 1"-y1mean)/y1sd
  Z_score_2 = (y$"LGR 2"-y1mean)/y1sd
  Z_score_3 = (y$"LGR 3"-y1mean)/y1sd
  Z_score_4 = (y$"LGR 4"-y1mean)/y1sd
  
  yz = mutate(y,Z_score_1,Z_score_2,Z_score_3,Z_score_4)
  
  ####################
  ####calculate mean and standard deviation of replictaes####
  if(smooth){
  aggregatedmean <-  aggregate(xLGR3 [,5:16], by = list(xLGR3$Plate, xLGR3$meanRow, xLGR3$meanColumn), FUN=mean, na.rm=TRUE) 
  }else{
  aggregatedmean <-  aggregate(xLGR3 [,5:15], by = list(xLGR3$Plate, xLGR3$meanRow, xLGR3$meanColumn), FUN=mean, na.rm=TRUE)   
  }
  if(smooth){
  aggregatedmean = rename(aggregatedmean, mean_usLGR = usLGR, mean_LGR = LGR, Plate = Group.1, Row = Group.2, Column = Group.3)
  }else{
  aggregatedmean = rename(aggregatedmean, mean_LGR = LGR, Plate = Group.1, Row = Group.2, Column = Group.3)
  }
  
  if(smooth){
  aggregatedsd <- aggregate(xLGR3 [,5:16], by = list(xLGR3$Plate, xLGR3$meanRow, xLGR3$meanColumn), FUN=sd, na.rm=TRUE)
  }else{
  aggregatedsd <- aggregate(xLGR3 [,5:15], by = list(xLGR3$Plate, xLGR3$meanRow, xLGR3$meanColumn), FUN=sd, na.rm=TRUE)  
  }
  if(smooth){
  aggregatedsd = rename(aggregatedsd, sd_usLGR = usLGR, sd_LGR = LGR, Plate = Group.1, Row = Group.2, Column = Group.3)
  }else{
  aggregatedsd = rename(aggregatedsd, sd_LGR = LGR, Plate = Group.1, Row = Group.2, Column = Group.3)
  }
  
  if(smooth){
  aggregatedsd <- aggregatedsd [c(1:3,12:13)]
  }else{
  aggregatedsd <- aggregatedsd [c(1:3,12)]
  }
  
  aggregated <- inner_join(aggregatedmean, aggregatedsd, by = c("Plate", "Row", "Column"))
  
  if(smooth){
  aggregated <- aggregated [c(1:13,16:17)]
  }else{
  aggregated <- aggregated [c(1:12,15)]
  }
  aggregated1 <- transform(aggregated, as.numeric(as.character(Plate)))
  
  ####### calculate Z-score over mean of replicates
  mean_Z_score = c()
  mean_Z_score = as.vector(scale(aggregated1$mean_LGR))
  aggregated2 = mutate(aggregated1, mean_Z_score)
  
  ########## if you ever reading this: you're almost done! so have a glass of prosecco and celebrate! rn it's friday, half past 5 and i want to go home...
  keyfile = read_delim("keyfile.txt", delim="\t")
  system(paste("perl", "keyfile_letters_to_numbers.pl"))
  
  
  newkeyfile = read_delim("newkeyfile.txt", delim="\t") 
  newkeyfile = rename(newkeyfile, Plate = `Plate #`)
  newkeyfile = transform(newkeyfile, Plate = as.numeric(Plate))
  
  mean_file <- inner_join(aggregated2, newkeyfile, by = c("Plate", "Row", "Column"))
  
  if(smooth){
    mean_file = mean_file[c(1:5,8:17)]
  }else{
    mean_file = mean_file[c(1:5,8:15)]
  }
  
  if(smooth){
    write.csv(mean_file[, c(1:3,15,4:14)],  file = sprintf("mean_file_%s.csv", Var))
  }else{
    write.csv(mean_file[, c(1:3,13,4:12)],  file = sprintf("mean_file_%s.csv", Var))
  }
  
  if(smooth){
  write.csv(yz[, c(1,3,4,2,5:8,29:32,21:28,13:20,9:12,33:40)],  file = sprintf("replicates_file_%s.csv", Var))
  }else{
  write.csv(yz[, c(1,3,4,2,5:8,29:32,21:28,13:20,9:12,33:36)],  file = sprintf("replicates_file_%s.csv", Var))
  }
  
  return(xrearranged)
}
