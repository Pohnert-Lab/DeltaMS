# plots isotopologue distributions from a label report to PDF named outputfile; 12 plots/page
# Highligts the plots of the selected isotopic ratio
# Distance between Peaks in isotopologue group (28.01.2016)
plotLabelReport_isoRat12 <- function(isoLabelReport, intOption = "rel",iRatio, errRatio, disPeak, maxNumDelta, classes, labeledSamples, outputfile) {
  # distance between peaks # 28.01.2016
  # new Parameter "maxNumDelta" and "disPeak" are missing in old implementation
  # Standard value setting: 13.07.2016
  if (missing(maxNumDelta)) {  # Number of Labels
    maxNumDelta <- 20;
  } 
  if (missing(disPeak))  { # IMD
    disPeak <- 10000;
  } 
  # isotopic ratio # 18.01.2016
  if (is.null(iRatio))  {
    iRatio <- c(0.5, 0.5);
  } 
  if (is.null(errRatio))  {
    errRatio <-0.05;
  } 
  # fixed configuration paramter
  # Color of the "correct" peaks   
  col1 <- "darkblue"
  col2 <- "blue"
  #  
  
  
  B1 <- iRatio[1]/(iRatio[1]+iRatio[2])
  B2 <- iRatio[2]/(iRatio[1]+iRatio[2])
  B_I1 <- c(0,0)
  B_I2 <- c(0,0)
  B_I1[1] <- B1*(1-errRatio)
  B_I1[2] <- B1*(1+errRatio)
  B_I2[1] <- B2*(1-errRatio)
  B_I2[2] <- B2*(1+errRatio)
  rat12 <- c(0,0)
  #  
  pdf(outputfile, width = 8.5, height = 11);
  numGroups = length(isoLabelReport[[1]]);
  par(mfrow=c(4,3));
  iz <- 1
  for (i in 1:numGroups) {
    if (iz %% 12 == 1) {
      par(mfrow=c(4,3));
    }
    if (intOption == "rel") {
      toPlot = cbind(isoLabelReport$meanRelU[[i]], isoLabelReport$meanRelL[[i]]);
      sds = cbind(isoLabelReport$sdRelU[[i]], isoLabelReport$sdRelL[[i]]);
      YLIM = c(0,1);
      rat12[1] <- toPlot[1,1]/(toPlot[1,1]+toPlot[2,2])
      rat12[2] <- toPlot[2,2]/(toPlot[1,1]+toPlot[2,2])
      ys <- 0.2
    }
    else if (intOption == "abs") {
      toPlot = cbind(isoLabelReport$meanAbsU[[i]], isoLabelReport$meanAbsL[[i]]);
      sds = cbind(isoLabelReport$cvTotalU[[i]]*isoLabelReport$meanAbsU[[i]], isoLabelReport$cvTotalL[[i]]*isoLabelReport$meanAbsL[[i]]);
      YLIM = NULL;
      totalPoolsU = colSums(isoLabelReport$sampleData[[i]][,which(classes != labeledSamples)]);
      totalPoolsL = colSums(isoLabelReport$sampleData[[i]][,which(classes == labeledSamples)]);
      T = t.test(totalPoolsU, totalPoolsL);
      pval = round(T$p.value, 3);
      rat12[1] <- toPlot[1,1]/(toPlot[1,1]+toPlot[2,2])
      rat12[2] <- toPlot[2,2]/(toPlot[1,1]+toPlot[2,2])
      ys <- max(toPlot) * 0.2
    }
    else {
      print("Error. Specify intOption as rel or abs.");
      break;
    }
    numPeaks = dim(toPlot)[1];
    
    diff21 <- isoLabelReport$isotopologue[[i]][2]- isoLabelReport$isotopologue[[i]][1]
    isDiff <- round(diff21)
    if (isDiff > round(disPeak*maxNumDelta)) {next};
    
    iz <- iz + 1
    rownames(toPlot) = as.character(round(isoLabelReport$isotopologue[[i]], 4));
    colnames(toPlot) = c("U", "L");
    
    if ((intOption == "rel") || (intOption == "abs")) {
      if (rat12[1] >= B_I1[1] && rat12[1] <= B_I1[2] && rat12[2] >= B_I2[1] && rat12[2] <= B_I2[2]) {
        mp = barplot(t(toPlot), beside = TRUE, legend.text = c("I1", "I2"), 
                     main = paste("m/z = ", rownames(toPlot)[1], "\n", "RT = ", isoLabelReport$rt[[i]][1]), 
                     axisnames = FALSE, ylim = YLIM , col = c(col1, col2)); #, col = c("dark red", "grey")
    #    text(seq(1.5, (numPeaks-1)*2+1.5, by = 2), pos=2, par("usr")[3]-0.2, xpd = TRUE, labels = rownames(toPlot), srt = 45);
              text(seq(1.5, (numPeaks-1)*3+1.5, by = 3), par("usr")[3]-ys, xpd = TRUE, labels = rownames(toPlot), srt = 45);
      } else {
        mp = barplot(t(toPlot), beside = TRUE, legend.text = c("I1", "I2"), 
                     main = paste("m/z = ", rownames(toPlot)[1], "\n", "RT = ", isoLabelReport$rt[[i]][1]), 
                     axisnames = FALSE, ylim = YLIM , col = c("dimgray", "gray")); #, col = c("dark red", "grey")
  #      text(seq(1.5, (numPeaks-1)*2+1.5, by = 2), pos=2, par("usr")[3]-0.2, xpd = TRUE, labels = rownames(toPlot), srt = 45);
        text(seq(1.5, (numPeaks-1)*3+1.5, by = 3), par("usr")[3]-ys, xpd = TRUE, labels = rownames(toPlot), srt = 45);
      }
      
#      mp = barplot(t(toPlot), beside = TRUE, legend.text = c("U", "L"), 
#                   main = paste("m/z = ", rownames(toPlot)[1], "\n", "RT = ", isoLabelReport$rt[[i]][1]), 
#                   axisnames = FALSE, ylim = YLIM);
#      text(seq(1.5, (numPeaks-1)*3+1.5, by = 3), par("usr")[3]-0.2, xpd = TRUE, labels = rownames(toPlot), srt = 45);
      
      segments(mp, t(toPlot) - t(sds), mp, t(toPlot) + t(sds), lwd = 2 );
    }
    else {
      barplot(toPlot, legend.text = rownames(toPlot), 
              main = paste("m/z = ", rownames(toPlot)[1], "\n", "RT = ", isoLabelReport$rt[[i]][1], "\n", "p = ", pval));
    }
  }
  ################################################
  # the remaining groups of peaks
  
  
  ################################################
  dev.off();
}
