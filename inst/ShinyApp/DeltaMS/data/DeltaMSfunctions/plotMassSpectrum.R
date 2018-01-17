# Jan 2017
# Plot MassSepectrum  
# plots isotopologue peaks from a label report to PDF named outputfile
# Highligts the plots of the selected isotopic ratio (?)
# Distance between Peaks in isotopologue group (28.01.2016)
plotMassSpectrum <- function(isoLabelReport, rawFiles, iRatio, errRatio, disPeak, maxNumDelta, RTwin, maxLT, infoList, classes, labeledSamples, labIso, outputfile) {
  # isoLabelReport = Output of getIsoLabelReport_doubleLabel().
  # rawFiles = Names of the raw files
  # iRatio = Ratio of the added isotopes
  # errRatio = Relative error of the isotope ratio
  # #  If the found peaks correspond to the predetermined isotope ratio, they are marked pink.
  # disPeak = distance between peaks # 28.01.2016 # Peaks with distance greater "disPeak*maxNumDelta" are not plotted
  # maxNumDelta = number of peaks 
  # RTwin = timeframe of RT 
  # #maxLT = Maximum number of peaks of a group to be plotted. # 14-02-2017 
  # classes = names of classes (intern)
  # labelSamples = class of the labeled sample (intern)
  if (is.null(maxNumDelta))  {
    maxNumDelta <- 10000;
  } 
  if (is.null(disPeak))  {
    disPeak <- 1;
  } 
  #  Maximum number of peaks to plot:  maxLT (14-02-17) 
  if (is.null(maxLT))  {
    maxLT <- 2;
  } 
  
  # isotopic ratio # 18.01.2016
  if (is.null(iRatio))  {
    iRatio <- c(0.5, 0.5);
  } 
  if (is.null(errRatio))  {
    errRatio <-0.05;
  } 
  if (is.null(RTwin))  {
    RTwin <-10;
  } 
  # legende # 10.05.2017
  if (is.null(labIso))  { 
    labIso <- c("I1", "I2");
  } 
  
  if (is.null(infoList))  {
    return("Error, Parameter 'infoList' is need! Perform first the function 'InfoOsoList'")
  } 
  
  className <-  as.matrix(classes)
  # fixed configuration paramter
  # Color of the "correct" peaks   
  col1 <- "darkblue"
  #  
  #  Classification of Isotopic pattern
  colI1 <- "darkorchid4"
  colI2 <- "darkorchid1"
  ## end: fixed ...  
  
  
  # Number of Isotopic pattern
  numPat <- length(iRatio)
  
  # Number of the calculated Isotopic pattern
  calcPat <- dim(infoList$iPat)[2]
  
  # Number useed
  curPat <- min(c(numPat, calcPat, maxLT))
  
  
  # Sort the data for the output
  outIDX <- order(infoList$infoMat[,1], infoList$infoMat[,2],infoList$infoMat[,3], infoList$infoMat[,4], decreasing = c(FALSE, FALSE, TRUE, TRUE), method="radix")
  
  # Count raw-files
  numRaw <- length(rawFiles)
  listRaw <- list()
  for (i in 1:numRaw){
    listRaw[[i]] <- xcmsRaw(rawFiles[i])
  }
  # cCol <- max(c(3, numRaw))
  cCol <- 3
  numPlot <- cCol * 5
  #
  pdf(outputfile, width = 8.5, height = 11);
  # x11()
  numGroups = length(isoLabelReport[[1]]);
 # par(mfrow=c(5,cCol));
#  iz <- 1
  

  # the remaining groups of peaks
  iz <- 1
  maxDist <- max(infoList$infoMat[,1]) # Maximum distance of the peaks (factor)
  
  for (ic in 1 : maxDist) {
    if (ic > curPat) {next}
    if (ic == 1) {
      #     par(mfrow=c(4,3));
      par(mfrow=c(5,cCol));
    }
    toPlot <- iRatio/sum(iRatio) #(as.matrix(iRatio[[1]]))
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

#  Format of the output-file new line or not ?     
#    for (k in 2:cCol) {
#      plot(c(0, 1), c(0, 1), type="l")
#      lines(c(0, 1), c(1, 0), type="l")
#      iz <- iz +1
#    }
    
    for (i in 1:numGroups) {
      # Each group starts with a new line
      if (!(iz %% cCol == 0)) {
        iz <- iz+ cCol - (iz %% cCol)
      } 
      if (iz %% numPlot == 1) {
        par(mfrow=c(5,cCol));
      }
      # current index
      cIdx <- outIDX[i]
      if (infoList$infoMat[cIdx,1] != ic) {
        next
      }
      # Determination of the data for the output
      toPlot = cbind(isoLabelReport$meanAbsU[[cIdx]]);
      sds = cbind(isoLabelReport$cvTotalU[[cIdx]]*isoLabelReport$meanAbsU[[cIdx]]);
      YLIM = c(0, max(toPlot));
      rat12 <- toPlot/sum(toPlot)
      
      #
      numPeaks = dim(toPlot)[1];
      
      mzPeaks <- isoLabelReport$isotopologue[[cIdx]]
      rtPeaks <- isoLabelReport$rt[[cIdx]]
      rownames(toPlot) = as.character(round(isoLabelReport$isotopologue[[cIdx]], 3));
      
      # RT-window as parameter ?    
      for (k in 1:numRaw) {
        iz <- iz + 1
        idRT <- which(listRaw[[k]]@scantime>rtPeaks[1]-0.05*RTwin & listRaw[[k]]@scantime<rtPeaks[numPeaks]+0.05*RTwin)
        if (length(idRT) < 1) {
          idRT <- which(listRaw[[k]]@scantime>rtPeaks[1]-0.2*RTwin & listRaw[[k]]@scantime<rtPeaks[numPeaks]+0.2*RTwin)
          if (length(idRT) < 1) {
            idRT <- which(listRaw[[k]]@scantime>rtPeaks[1]-0.5*RTwin & listRaw[[k]]@scantime<rtPeaks[numPeaks]+0.5*RTwin)
          }
          if (length(idRT) < 1) {
            idRT <- which(listRaw[[k]]@scantime>rtPeaks[1]-0.75*RTwin & listRaw[[k]]@scantime<rtPeaks[numPeaks]+0.75*RTwin)
          }
        }
        if (length(idRT) < 1) {
          plot(c(0, 1), c(0, 1), type="l")
          lines(c(0, 1), c(1, 0), type="l")
          next
        }
        df <- disPeak
        
        ## handle last spectrum
        if (idRT[1] == length(listRaw[[k]]@scanindex)) {
          followingScanIndex <- length(listRaw[[k]]@env$mz)
        } else {
          followingScanIndex <- listRaw[[k]]@scanindex[idRT[1]+1]
        }
        
        idx <- (listRaw[[k]]@scanindex[idRT[1]]+1):min(followingScanIndex,
                                                       length(listRaw[[k]]@env$mz), na.rm=TRUE)
        mzrange <- c(mzPeaks[1]-df, mzPeaks[numPeaks]+df)
        if (length(mzrange) >= 2) {
          mzrange <- range(mzrange)
          idx <- idx[listRaw[[k]]@env$mz[idx] >= mzrange[1] & listRaw[[k]]@env$mz[idx] <= mzrange[2]]
        }
        
        points <- cbind(listRaw[[k]]@env$mz[idx], listRaw[[k]]@env$intensity[idx])
        title = paste("#", i, ": ", className[k],"\n", "#",k, ": ", round(listRaw[[k]]@scantime[idRT[1]]/60, 2),
                      " min (scan ",idRT[1], ")", sep = "")
 #       plot(points, type="h", main = title, xlab="m/z", ylab="Intensity")
        if (length(idx) > 0) {
          plot(points, type="h", main = title, xlab="m/z", ylab="Intensity")
        } else {
          plot(c(0, 1), c(0, 1), type="l", main = title, xlab="m/z", ylab="Intensity")
          lines(c(0, 1), c(1, 0), type="l")
        }
        
        scanVal <- getScan(listRaw[[k]], idRT[1], mzrange = c(mzPeaks[1]-df, mzPeaks[numPeaks]+df))
        if (length(scanVal)==0){
          next
        }
        
#        hv1 <- scanVal[,1]- mzPeaks[1]
#        i1 <- scanVal[which.min(abs(hv1)),2]
#        hv2 <- scanVal[,1]- mzPeaks[2]
#        i2 <- scanVal[which.min(abs(hv2)),2]
        iA <- matrix(data=0, nrow=1, ncol=numPeaks)
        for (k in 1:numPeaks){
          hvA <- scanVal[,1]- mzPeaks[k]
          iA[k] <- scanVal[which.min(abs(hvA)),2]
          
        }
###        for (k in 1: numPeaks) {
###          text(mzPeaks[k], iA[k], rownames(toPlot)[k])
###          lines(c(mzPeaks[k], mzPeaks[k]), c(0,iA[k]),  type="l", lwd = 1)
          if (infoList$infoMat[cIdx,3] == 1) {
            for (k in 1: numPeaks) {
              text(mzPeaks[k], iA[k], rownames(toPlot)[k], col = c(col1))
              lines(c(mzPeaks[k], mzPeaks[k]), c(0,iA[k]),  col = c(col1), type="l", lwd = 2)
            }
          } else {
            for (k in 1: numPeaks) {
              text(mzPeaks[k], iA[k], rownames(toPlot)[k])
              lines(c(mzPeaks[k], mzPeaks[k]), c(0,iA[k]),  type="l", lwd = 1)
            }
          }
          
          
          
        
      }     
    }
  }
  
  
  
  ################################################
  
  dev.off();
}
