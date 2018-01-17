# plots isotopologue distributions from a label report to PDF named outputfile; 12 plots/page
# Highligts the plots of the selected isotopic ratio
# Distance between Peaks in isotopologue group (28.01.2016)
plotLabelReport_isoRat <- function(isoLabelReport, intOption = "rel",iRatio, errRatio, disPeak, maxNumDelta, maxLT, infoList,  classes, labeledSamples, labIso, outputfile) {
  # distance between peaks # 28.01.2016
  # new Parameter "maxNumDelta" and "disPeak" are missing in old implementation
  # Standard value setting: 13.07.2016
  if (missing(maxNumDelta)) {  # maximum number of DELTA's 
    maxNumDelta <- 20;
  } 
  if (missing(disPeak))  { # IMD, distance between two peaks
    disPeak <- 1;
  } 
  # isotopic ratio # 18.01.2016
  if (is.null(iRatio))  {
    iRatio <- c(0.5, 0.5);
  } 
  if (is.null(errRatio))  {
    errRatio <-0.05;
  } 
  #  Maximum number of peaks to plot:  maxLT (14-02-17) 
  if (is.null(maxLT))  {
    maxLT <- 2;
  } 
  # legende # 10.05.2017
  if (is.null(labIso))  { 
    labIso <- c("I1", "I2");
  } 
  
  if (is.null(infoList))  {
    return("Error, Parameter 'infoList' is need! Perform first the function 'InfoOsoList'")
  } 
  
  # fixed configuration paramter
  # Color of the "correct" peaks   
  col1 <- "darkblue"
  col2 <- "blue"
  #  
  colN1 <- "dimgray"
  colN2 <- "gray"
  #  Classification of Isotopic pattern
  colI1 <- "darkorchid4"
  colI2 <- "darkorchid1"
  ## end: fixed ...  
  
  
 # B1 <- iRatio[1]/(iRatio[1]+iRatio[2])
#  B2 <- iRatio[2]/(iRatio[1]+iRatio[2])
#  B_I1 <- c(0,0)
#  B_I2 <- c(0,0)
#  B_I1[1] <- B1*(1-errRatio)
#  B_I1[2] <- B1*(1+errRatio)
#  B_I2[1] <- B2*(1-errRatio)
#  B_I2[2] <- B2*(1+errRatio)
#  rat12 <- c(0,0)
  #  
  iRatio <- iRatio/sum(iRatio)
  # 
  numGroups = length(isoLabelReport[[1]]); # Number of ISO - groups
  # Number of isotopic pattern 
  numPat <- 1 # always one if samples labled with one isotope # length(iRatio)
  
  interMax <- 1
  # Intervall of interest for the isotopologue pattern 
  for (i in 1:numGroups) {
    df <- isoLabelReport$isotopologue[[i]]
    ldf <- length(df)
    fDist <- round((df[ldf] - df[1])/disPeak)
    if (fDist > interMax){
      interMax <- fDist
    }
  }
  
  curInt <- min(c(interMax, maxNumDelta))
  
  numPlot <- max(infoList$infoMat[,1]) # Maximum distance of the peaks (factor)
  
  # Sort the data for the output
  outIDX <- order(infoList$infoMat[,1], infoList$infoMat[,2],infoList$infoMat[,3], infoList$infoMat[,4], decreasing = c(FALSE, FALSE, TRUE, TRUE), method="radix")
  
  
  #
  pdf(outputfile, width = 8.5, height = 11);
  numGroups = length(isoLabelReport[[1]]);
  par(mfrow=c(4,3));
  iz <- 1
  for (ic in 1 : numPlot){
    if (ic > curInt) {
      next
    }
    toPlot <- iRatio
    rWord <-  c("m/z")
    rNames <- c("m/z")
    
    for (i in 1) {
      rNames <- cbind(rNames, paste(rWord, "+", format(disPeak*ic, digits = 4)) )
    }
    
    YLIM <- c(0, 1)
    mp = barplot(t(toPlot), beside = TRUE, legend.text = labIso, 
                 main = paste("Distance ", ic, "\n", rNames[2]), 
                 axisnames = FALSE, ylim = YLIM , col = c(colI1, colI2)); #, col = c("dark red", "grey")
    text(seq(1.5, (ic)*2+1.5, by = 2), pos=2, par("usr")[3]-0.15, xpd = TRUE, labels = rNames, srt = 45);
    iz <- iz + 1
#    if (ic == 5){
#      print(ic)
#    }
    
    for (i in 1:numGroups) {
      if (iz %% 12 == 1) {
          par(mfrow=c(4,3));
      }
      # current index
      cIdx <- outIDX[i]
      if (infoList$infoMat[cIdx,1] != ic) {
        next
      }
      
      if (intOption == "rel") {
          toPlot = cbind(isoLabelReport$meanRelU[[cIdx]], isoLabelReport$meanRelL[[cIdx]]);
          sds = cbind(isoLabelReport$sdRelU[[cIdx]], isoLabelReport$sdRelL[[cIdx]]);
          YLIM = c(0,1);
    #      rat12[1] <- toPlot[1,1]/(toPlot[1,1]+toPlot[2,2])
    #      rat12[2] <- toPlot[2,2]/(toPlot[1,1]+toPlot[2,2])
          ys <- 0.2
      }
      else if (intOption == "abs") {
         toPlot = cbind(isoLabelReport$meanAbsU[[cIdx]], isoLabelReport$meanAbsL[[cIdx]]);
         sds = cbind(isoLabelReport$cvTotalU[[cIdx]]*isoLabelReport$meanAbsU[[cIdx]], isoLabelReport$cvTotalL[[cIdx]]*isoLabelReport$meanAbsL[[cIdx]]);
         YLIM = NULL;
     #    totalPoolsU = colSums(isoLabelReport$sampleData[[cIdx]][,which(classes != labeledSamples)]);
    #     totalPoolsL = colSums(isoLabelReport$sampleData[[cIdx]][,which(classes == labeledSamples)]);
    #     T = t.test(totalPoolsU, totalPoolsL);
    #     pval = round(T$p.value, 3);
    #     rat12[1] <- toPlot[1,1]/(toPlot[1,1]+toPlot[2,2])
    #     rat12[2] <- toPlot[2,2]/(toPlot[1,1]+toPlot[2,2])
         ys <- max(toPlot) * 0.2
      }
      else {
         print("Error. Specify intOption as rel or abs.");
         break;
      }
      numPeaks = dim(toPlot)[1];
    
      diff21 <- isoLabelReport$isotopologue[[cIdx]][2]- isoLabelReport$isotopologue[[cIdx]][1]
      isDiff <- round(diff21)
      if (isDiff > round(disPeak*maxNumDelta)) {next};
    
      
      rownames(toPlot) = as.character(round(isoLabelReport$isotopologue[[cIdx]], 4));
      colnames(toPlot) = c("U", "L");
      
      aDim <- dim(toPlot)
      if (aDim[1] > maxLT){
        toPlot <-  toPlot[-c((maxLT+1):aDim[1]), ]
        sds <-  sds[-c((maxLT+1):aDim[1]), ]
      }
    
      if ((intOption == "rel") || (intOption == "abs")) {
          iz <- iz + 1
#         if (rat12[1] >= B_I1[1] && rat12[1] <= B_I1[2] && rat12[2] >= B_I2[1] && rat12[2] <= B_I2[2]) {
           if (infoList$infoMat[cIdx,3] == 1) {
             if (infoList$infoMat[cIdx,4] == 1) {
                mp = barplot(t(toPlot), beside = TRUE, legend.text = labIso, 
                     main = paste("m/z = ", rownames(toPlot)[1], "\n", "RT = ", round((isoLabelReport$rt[[cIdx]][1]/60),2)), 
                     axisnames = FALSE, ylim = YLIM , col = c(col1, col2)); #, col = c("dark red", "grey")
             } else {
               aCol <- matrix(data=NA, nrow = 1, ncol = (2*numPeaks))
               for (k in 1:numPeaks){
                 if (k<infoList$infoMat[cIdx,4]){
                   aCol[2*k-1] <- colN1
                   aCol[2*k] <- colN2
                 } else {
                   aCol[2*k-1] <- col1
                   aCol[2*k] <- col2
                 }
               }
                mp = barplot(t(toPlot), beside = TRUE, legend.text = labIso, 
                          main = paste("m/z = ", rownames(toPlot)[1], "\n", "RT = ", round((isoLabelReport$rt[[cIdx]][1]/60),2)), 
                          axisnames = FALSE, ylim = YLIM , col = aCol); #, col = c("dark red", "grey")
             }
             #    text(seq(1.5, (numPeaks-1)*2+1.5, by = 2), pos=2, par("usr")[3]-0.2, xpd = TRUE, labels = rownames(toPlot), srt = 45);
            text(seq(1.5, (numPeaks-1)*3+1.5, by = 3), par("usr")[3]-ys, xpd = TRUE, labels = rownames(toPlot), srt = 45);
           } else {
             mp = barplot(t(toPlot), beside = TRUE, legend.text = labIso, 
                     main = paste("m/z = ", rownames(toPlot)[1], "\n", "RT = ", round((isoLabelReport$rt[[cIdx]][1]/60),2)), 
                     axisnames = FALSE, ylim = YLIM , col = c("dimgray", "gray")); #, col = c("dark red", "grey")
        #        text(seq(1.5, (numPeaks-1)*2+1.5, by = 2), pos=2, par("usr")[3]-0.2, xpd = TRUE, labels = rownames(toPlot), srt = 45);
             text(seq(1.5, (numPeaks-1)*3+1.5, by = 3), par("usr")[3]-ys, xpd = TRUE, labels = rownames(toPlot), srt = 45);
           }
      
       #      mp = barplot(t(toPlot), beside = TRUE, legend.text = c("U", "L"), 
      #                   main = paste("m/z = ", rownames(toPlot)[1], "\n", "RT = ", round((isoLabelReport$rt[[cIdx]][1]/60),2)), 
      #                   axisnames = FALSE, ylim = YLIM);
      #      text(seq(1.5, (numPeaks-1)*3+1.5, by = 3), par("usr")[3]-0.2, xpd = TRUE, labels = rownames(toPlot), srt = 45);
      
         segments(mp, t(toPlot) - t(sds), mp, t(toPlot) + t(sds), lwd = 2 );
       }
       else {
         barplot(toPlot, legend.text = rownames(toPlot), 
               main = paste("m/z = ", rownames(toPlot)[1], "\n", "RT = ", round((isoLabelReport$rt[[cIdx]][1]/60),2), "\n", "p = ", pval));
       } 
     }
  }
  
  ################################################
  dev.off();
}
