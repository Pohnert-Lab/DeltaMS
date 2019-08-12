# Jan 2017
# Plot EIC  
# plots isotopologue peaks from a label report to PDF named outputfile
# Highligts the plots of the selected isotopic ratio (?)
# Distance between Peaks in isotopologue group (28.01.2016)
plotEIC_doubleLabel <- function(isoLabelReport, rawFiles, iRatio, errRatio, disPeak, maxNumDelta, RTwin, maxLT, infoList,  classes, labIso, labeledSamples, outputfile) {
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
  
  className <-  as.matrix(classes)
  
  # fixed configuration paramter
  # Color of the "correct" peaks   
  col1 <- "darkblue"
  col0 <- "black"
  #  Color: Classification of Isotopic pattern
  colI1 <- "darkorchid4"
  colI2 <- "darkorchid1"
  # Color of the "not-correct" peaks   
  colN1 <- "gray40"
  colN2 <- "gray80" # "lightblue" #
  
  
  # Color of the "correct" peaks   
  col1 <- "darkblue"
  col2 <- "blue" # "lightblue" #
  #  Classification of Isotopic pattern
  colI1 <- "darkorchid4"
  colI2 <-  "orchid1" #"darkorchid1"
  
  colfuncN <- colorRampPalette(c(colN1, colN2))
  colfuncI <- colorRampPalette(c(colI1, colI2))
  colfunc <- colorRampPalette(c(col1, col2))
  
  
  #
  dMZ <- 0.01   # aboslutes m/z-Intervall
  dRT <- 5      # factor for RT window
  #  end fixed
  
  # Number of Isotopic pattern
  numPat <- length(iRatio)
  
  # Number of the calculated Isotopic pattern
  calcPat <- dim(infoList$iPat)[2]
  
  # Number useed
  curPat <- min(c(numPat, calcPat, maxLT))
  
  # Sort the data for the output
  outIDX <- order(infoList$infoMat[,1], infoList$infoMat[,2],infoList$infoMat[,3], infoList$infoMat[,4],infoList$infoMat[,5], decreasing = c(FALSE, FALSE, TRUE, TRUE,TRUE), method="radix")
  
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
  par(mfrow=c(5,cCol));
  iz <- 1
  # for all isotopic patterns
  for (ic in 1 : numPat) {
    if (ic > curPat) { #Currently possible number
      break
    }
    if (ic == 1) {
      #     par(mfrow=c(4,3));
      par(mfrow=c(5,cCol));
    }
    toPlot <- (as.matrix(iRatio[[ic]]))
    rWord <-  c("m/z")
    rNames <- c("m/z")
    legTxt <- paste(ic,"x",labIso[1])
    
    
    legPat <- vector()
    for (i in 1:(ic)) {
      rNames <- cbind(rNames, paste(rWord, "+", format(disPeak*i, digits = 4)) )
      legPat[i] <- paste0("P", i)
      if (i == ic) {
        legTxt <- cbind(legTxt, paste(ic,"x",labIso[2]))
      }else{
        legTxt <- cbind(legTxt, paste(ic-i,"x",labIso[1], ",", i,"x",labIso[2]))
      }
    }
    legPat[ic+1] <- paste0("P", ic+1)
    
    # Output of the pattern to be identified
    YLIM <- c(0, 1)
    mp = barplot(t(toPlot), beside = TRUE, legend.text = legTxt,
                 main = paste("Level ", ic, "\n", "Isotopic pattern"), 
                 axisnames = FALSE, ylim = YLIM , col = colfuncI(ic+1));  #, col = c("dark red", "grey")
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
      

      #
      # Determination of the data for the output
      toPlot = cbind(isoLabelReport$meanAbsU[[cIdx]]);
      sds = cbind(isoLabelReport$cvTotalU[[cIdx]]*isoLabelReport$meanAbsU[[cIdx]]);
      YLIM = c(0, max(toPlot));
      rat12 <- toPlot/sum(toPlot)
      
      numPeaks = dim(toPlot)[1];

      mzPeaks <- isoLabelReport$isotopologue[[cIdx]]
      rtPeaks <- isoLabelReport$rt[[cIdx]]
      rownames(toPlot) = as.character(round(isoLabelReport$isotopologue[[cIdx]], 3));
    
    # RT-window as parameter ?    
      for (k in 1:numRaw) {
        if ((((infoList$infoMat[cIdx,1] == 1) && infoList$infoMat[cIdx,2] == (ic+1) && infoList$infoMat[cIdx,4] == 1))
            || (((infoList$infoMat[cIdx,1] == 1) && infoList$infoMat[cIdx,5] == (ic+1) ))) {
          
            if ((infoList$infoMat[cIdx,3] == 1)) { #infoList$iPat[cIdx,ic] == 1)) {
               colA <- col1
            } else {
               colA <- col0
            }
            pointsAll <- list()
            maxInt <- matrix(data=0, nrow = 1, ncol=numPeaks)
            for (m in 1:numPeaks){
              EIC1 <-  rawEIC(listRaw[[k]], mzrange = c(mzPeaks[m]- dMZ ,mzPeaks[m]+ dMZ), rtrange = c(rtPeaks[1]-dRT*RTwin, rtPeaks[numPeaks]+dRT*RTwin))
              maxInt[m] <- max(EIC1$intensity)
              pointsAll[[m]]<- cbind((listRaw[[k]]@scantime[EIC1$scan])/60, EIC1$intensity)
              
            }
            Imax <- max(maxInt)
            titleE <- paste("#", i, ": ", className[k],"\n", "#", k,  " EIC m/z  ",round(mzPeaks[1], 4)," / ", round(mzPeaks[numPeaks], 4),sep="")
            plot(pointsAll[[1]], type="l", main=titleE, xlab="min",
                        ylab="Intensity", xlim=c((rtPeaks[1]-5*RTwin)/60, (rtPeaks[2]+5*RTwin)/60), ylim = c(0, Imax), col = colA)
            for (m in 2:numPeaks) {
                lt <- m
                if (m > 6) {
                  lt <- 6
                }
                points(pointsAll[[m]], type ="l", lty = lt, col = colA)
            }
            if (ic == 1) {
              legend("topright", legend = labIso, col = colA, lty = 1:2, merge = TRUE)   
            }  else {
              legend("topright", legend = legPat, col = colA, lty = 1:(ic+1), merge = TRUE)   
            }
        }

      }
    
    }
  }
  
  
  ######################################
  # remaining groups of peaks
  notUse <- 1 # makes no sense
  if (!notUse) { 
  iz <- 1
  maxDist <- max(infoList$infoMat[,1]) # Maximum distance of the peaks (factor)
  
  for (ic in 2 : maxDist) {
    if (ic == 2) {
      #     par(mfrow=c(4,3));
      par(mfrow=c(5,cCol));
    }
    toPlot <- (as.matrix(iRatio[[1]]))
    rWord <-  c("m/z")
    rNames <- c("m/z")
    
    for (i in 1) {
      rNames <- cbind(rNames, paste(rWord, "+", format(disPeak*ic, digits = 4)) )
    }
    
    YLIM <- c(0, 1)
    mp = barplot(t(toPlot), beside = TRUE, legend.text = c("I1", "I2"), 
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
      
      numPeaks = dim(toPlot)[1];

      mzPeaks <- isoLabelReport$isotopologue[[cIdx]]
      rtPeaks <- isoLabelReport$rt[[cIdx]]
      rownames(toPlot) = as.character(round(isoLabelReport$isotopologue[[cIdx]], 3));
      
      # RT-window as parameter ?   
      colA <- col0
      for (k in 1:numRaw) {

        pointsAll <- list()
        maxInt <- matrix(data=0, nrow = 1, ncol=numPeaks)
        for (m in 1:numPeaks){
          EIC1 <-  rawEIC(listRaw[[k]], mzrange = c(mzPeaks[m]- dMZ ,mzPeaks[m]+ dMZ), rtrange = c(rtPeaks[1]-dRT*RTwin, rtPeaks[numPeaks]+dRT*RTwin))
          maxInt[m] <- max(EIC1$intensity)
          pointsAll[[m]]<- cbind((listRaw[[k]]@scantime[EIC1$scan])/60, EIC1$intensity)
          
        }
        Imax <- max(maxInt)
        title <- paste("#", i, ": ", className[k],"\n", "#", k,  " EIC m/z  ",round(mzPeaks[1], 1)," - ", round(mzPeaks[numPeaks], 1),sep="")
        plot(pointsAll[[1]], type="l", main=title, xlab="min",
             ylab="Intensity", xlim=c((rtPeaks[1]-5*RTwin)/60, (rtPeaks[2]+5*RTwin)/60), ylim = c(0, Imax), col = colA)
        for (m in 2:numPeaks) {
          lt <- m
          if (m > 6) {
            lt <- 6
          }
          points(pointsAll[[m]], type ="l", lty = lt, col = colA)
        }
        
        
      }
      
    }
  }
  } # notUse
  
  dev.off();
}
