IsoRatioCalc<-function(iRatioPerc, m = 2, mzdiff = 1){
  if(m < 2){
    return("Error, value has to be at least 2")
  }
  ## Get Matrix with abundances, consists of all isotopic combinations regarding m with replacement, but unordered.
  abundCombReplace<-CombSet(iRatioPerc,m,repl = TRUE, ord = FALSE) 
  if (iRatioPerc[1] == iRatioPerc[2]){
    tMat <- matrix(data=0, nrow = (m+1), ncol=m)
    for (i in 1:(m+1)){
      tMat[i,]<-abundCombReplace[1,]
    }
    abundCombReplace <- tMat
  }
    
  ## get matrix with all possible combinations with replacement
  freqVarReplace<-CombSet(c(0,1),m,repl = TRUE, ord = TRUE) 
  ## calculate deltaM values for all combinations, rowwise
  freqVarReplace2<-apply(freqVarReplace,1,sum) 
  #get table out this to show frequencies of calculated deltaM values 
  dfFreqVarReplace<-as.data.frame(table(freqVarReplace2))
  #Combine frequencies and abundancies matrix 
  abundCombReplace2<-cbind(abundCombReplace, dfFreqVarReplace$Freq)
  # multiply rowwise to get abundance for all deltaD values
  abundCombReplace3<-apply(abundCombReplace2,1,prod)
  #abundCombReplace4<-abundCombReplace3/max(abundCombReplace3)
  return(abundCombReplace3/sum(abundCombReplace3))
  # For illustration of ratios, these lines will generate a plot with isotopic patterns
 ## plot(abundCombReplace4,ylab = "Abundance [%]",xlab = "m/z + X", type = "h", lwd = 3, axes = FALSE, ylim = c(0,1))
##  axis(side = 1,at = seq(0:length(abundCombReplace4)), labels = seq(0,length(abundCombReplace4)*mzdiff, by = mzdiff))
##  axis(side = 2, at = seq(0,1 , by = 0.2), labels = seq(0,1, by = 0.2))
}
