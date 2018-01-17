# plots isotopologue distributions from a label report to PDF named outputfile; 12 plots/page
# Highligts the plots of the selected isotopic ratio
# Distance between Peaks in isotopologue group (28.01.2016)
# plotLabelReport <- function(isoLabelReport, intOption = "rel", classes, labeledSamples, outputfile)
plotLabelReport_doubleLabel <- function(isoLabelReport, intOption = "rel", iRatio, errRatio, disPeak, maxNumDelta, maxLT = 2, infoList ,classes, labeledSamples, labIso, outputfile) {
    # Changes from the original X13CMS
  # iRatio = List of Isotopic pattern with their propabilities # 30-03-2017
  # errRatio = Relative error of the isotope ratio
  # #  If the found peaks correspond to the predetermined isotope ratio, they are marked pink.
  # disPeak = distance between peaks # 28.01.2016 # Peaks with distance greater "disPeak*maxNumDelta" are not plotted
  # maxLT = Maximum number of peaks of a group to be plotted. # 14-02-2017 
  # infoList = Summary of the procedure of pattern detection # 04-04-2017
  # labIso = Labels of Isotopes; used as legend
  #           iDist = matrix(number of found groups, Number of peaks in the distance of the mass difference of the isotopes )
  #                   1 = peak, NA = no peak
  #           iPat  = matrix similiar to iDist, Number of columns reduced by one
  #                   1 = The neighbouring peaks have the expected intensity ratio; 0 = not
  #           infoMat:  1. column = Distance between the first and second peak
  #                     2. column = Number of peaks in the found pattern
  #                     3. column =  Pattern meets all conditions -> 1
  #                     4. column = The greatest distance between adjacent peaks corresponds to the "simple" distance -> 1
  if (is.null(maxNumDelta))  {
    maxNumDelta <- 10000;
  } 
  #  Maximum number of peaks to plot:  maxLT (14-02-17) 
  if (is.null(maxLT))  {
    maxLT <- 2;
  } 
 
  if (is.null(disPeak))  {
    disPeak <- 1;
  } 
  
  # isotopic ratio # 18.01.2016
  if (is.null(iRatio))  {
    iRatio[[1]] <- c(0.5, 0.5);
  } 
  if (is.null(errRatio))  {
    errRatio <-0.05;
  } 
  # legende # 10.05.2017
  if (is.null(labIso))  { 
    labIso <- c("I1", "I2");
  } 
  if (is.null(infoList))  {
    return("Error, Parameter 'infoList' is need! Perform first the function 'InfoOsoList'")
  } 
  # Fix configurations parameter # 07-02-2017
  bMin = TRUE   # RT in min or sec
  if (bMin) {
    tFac = 60;  # RT in minutes
  }  else {
    tFac = 1;   # RT in seconds
  }
  
  # Color of the "not-correct" peaks   
  colN1 <- "gray40"
  colN2 <- "gray80" # "lightblue" #
  
   
  # Color of the "correct" peaks   
  col1 <- "darkblue"
  col2 <- "blue" # "lightblue" #
  #  Classification of Isotopic pattern
  colI1 <- "darkorchid4"
  colI2 <-  "orchid1" #"darkorchid1"
  
  bLeg <- 0 # legend Y/N
  ## end fix configuration
  
  colfuncN <- colorRampPalette(c(colN1, colN2))
  colfuncI <- colorRampPalette(c(colI1, colI2))
  colfunc <- colorRampPalette(c(col1, col2))
    #  
  
  # Number of Isotopic pattern
  numPat <- length(iRatio)
  
  # Number of the calculated Isotopic pattern
  calcPat <- dim(infoList$iPat)[2]
  
  # Number useed
  curPat <- min(c(numPat, calcPat, maxLT))

  
 # Sort the data for the output
  outIDX <- order(infoList$infoMat[,1], infoList$infoMat[,2],infoList$infoMat[,3], infoList$infoMat[,4],infoList$infoMat[,5], decreasing = c(FALSE, FALSE, TRUE, TRUE,TRUE), method="radix")
  
  #  
  pdf(outputfile, width = 8.5, height = 11);
  numGroups = length(isoLabelReport[[1]]);
  iz <- 1
  
  # for all isotopic patterns
  for (ic in 1 : numPat) {
    if (ic > curPat) { #Currently possible number
      break
    }
    if (ic == 1) {
       par(mfrow=c(4,3));
    }
    toPlot <- (as.matrix(iRatio[[ic]]))
    rWord <-  c("m/z")
    rNames <- c("m/z")
    legTxt <- paste(ic,"x",labIso[1])
    
    for (i in 1:(ic)) {
     rNames <- cbind(rNames, paste(rWord, "+", format(disPeak*i, digits = 4)) )
     if (i == ic) {
         legTxt <- cbind(legTxt, paste(ic,"x",labIso[2]))
     }else{
       legTxt <- cbind(legTxt, paste(ic-i,"x",labIso[1], ",", i,"x",labIso[2]))
     }
    }
   
    # Output of the pattern to be identified
    YLIM <- c(0, 1)
    mp = barplot(t(toPlot), beside = TRUE, legend.text = legTxt, 
                 main = paste("Level ", ic, "\n", "Isotopic pattern"), 
                 axisnames = FALSE, ylim = YLIM , col =  colfuncI(ic+1)); #c(colI1, colI2)); #, col = c("dark red", "grey")
    text(seq(1.5, (ic)*2+1.5, by = 2), pos=2, par("usr")[3]-0.15, xpd = TRUE, labels = rNames, srt = 45);
    iz <- iz + 1
    
    for (i in 1:numGroups) {
      if (iz %% 12 == 1) {
         par(mfrow=c(4,3));
      }
      cIdx <- outIDX[i]
      
    # Determination of the data for the output
    if (intOption == "rel") {
      toPlot = cbind(isoLabelReport$meanRelU[[cIdx]]);
      sds = cbind(isoLabelReport$sdRelU[[cIdx]]);
      YLIM = c(0,1);
      rat12 <- toPlot/sum(toPlot)
    }
    else if (intOption == "abs") {
      toPlot = cbind(isoLabelReport$meanAbsU[[cIdx]]);
      sds = cbind(isoLabelReport$cvTotalU[[i]]*isoLabelReport$meanAbsU[[cIdx]]);
      YLIM = c(0, max(toPlot));
      rat12 <- toPlot/sum(toPlot)
    }
    else {
      print("Error. Specify intOption as rel or abs.");
      break;
    }
   numPeaks = dim(toPlot)[1];

    
  
    rownames(toPlot) = as.character(round(isoLabelReport$isotopologue[[cIdx]], 4));
#   Change 2017: Variable number of peaks.
    lt <- length(toPlot)

    if ((intOption == "rel") || (intOption == "abs") ) {
      if ((((infoList$infoMat[cIdx,1] == 1) && infoList$infoMat[cIdx,2] == (ic+1) && infoList$infoMat[cIdx,4] == 1))
          || (((infoList$infoMat[cIdx,1] == 1) && infoList$infoMat[cIdx,5] == (ic+1) ))) {
         iz <- iz + 1
 #        if ((infoList$infoMat[cIdx,3] == 1) || (infoList$iPat[cIdx,ic] == 1)) {
         if ((infoList$infoMat[cIdx,3] >= 1) ) { # The search is not performed according to the current pattern in larger groups. 
             mp = barplot(t(toPlot), beside = TRUE, # legend.text = labIso, 
                   main = paste("m/z = ", rownames(toPlot)[1], "\n", "RT = ", format(isoLabelReport$rt[[cIdx]][1]/tFac, digits = 5)), 
                   axisnames = FALSE, ylim = YLIM , col = colfunc(lt)); #c(col1, col2)); #, col = c("dark red", "grey")
            text(seq(1.5, (numPeaks-1)*2+1.5, by = 2), pos=2, par("usr")[3]-0.2, xpd = TRUE, labels = rownames(toPlot), srt = 45);
        } else {
            mp = barplot(t(toPlot), beside = TRUE, # legend.text = labIso, 
                   main = paste("m/z = ", rownames(toPlot)[1], "\n", "RT = ", format(isoLabelReport$rt[[cIdx]][1]/tFac, digits = 5)), 
                   axisnames = FALSE, ylim = YLIM , col = colfuncN(lt)); #c("dimgray", "gray")); #, col = c("dark red", "grey")
            text(seq(1.5, (numPeaks-1)*2+1.5, by = 2), pos=2, par("usr")[3]-0.2, xpd = TRUE, labels = rownames(toPlot), srt = 45);
      
        }
        segments(mp, t(toPlot) - t(sds), mp, t(toPlot) + t(sds), lwd = 2 );
       }
    } 
    else {
      barplot(toPlot, legend.text = rownames(toPlot), 
              main = paste("m/z = ", rownames(toPlot)[1], "\n", "RT = ", isoLabelReport$rt[[i]][1], "\n", "p = ", pval));
    }
  }
  }
  ################################################
  # the remaining groups of peaks
  notUse <- 1
  if (!notUse) {
  iz <- 1
  numPlot <- max(infoList$infoMat[,1]) # Maximum distance of the peaks (factor)
  
  for (ic in 2 : numPlot) {
    if (ic == 2) {
      par(mfrow=c(4,3));
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
    
    for (i in 1:numGroups) {
      if (iz %% 12 == 1) {
        par(mfrow=c(4,3));
      }
      cIdx <- outIDX[i]
      if (infoList$infoMat[cIdx,1] != ic) {
        next
      }
        
     
      # Determination of the data to be showed
      if (intOption == "rel") {
        toPlot = cbind(isoLabelReport$meanRelU[[cIdx]]);
        sds = cbind(isoLabelReport$sdRelU[[cIdx]]);
        YLIM = c(0,1);
        rat12 <- toPlot/sum(toPlot)
      }
      else if (intOption == "abs") {
        toPlot = cbind(isoLabelReport$meanAbsU[[cIdx]]);
        sds = cbind(isoLabelReport$cvTotalU[[i]]*isoLabelReport$meanAbsU[[cIdx]]);
        YLIM = c(0, max(toPlot));
         rat12 <- toPlot/sum(toPlot)
      }
      else {
        print("Error. Specify intOption as rel or abs.");
        break;
      }
      numPeaks = dim(toPlot)[1];

      rownames(toPlot) = as.character(round(isoLabelReport$isotopologue[[cIdx]], 4));
      #   Change 2017: Variable number of peaks.
      lt <- length(toPlot)

     if ((intOption == "rel") || (intOption == "abs") ) {
         mp = barplot(t(toPlot), beside = TRUE, legend.text = c("I1", "I2"), 
                      main = paste("m/z = ", rownames(toPlot)[1], "\n", "RT = ", format(isoLabelReport$rt[[cIdx]][1]/tFac, digits = 5)), 
                       axisnames = FALSE, ylim = YLIM , col = c("dimgray", "gray")); #, col = c("dark red", "grey")
         text(seq(1.5, (numPeaks-1)*2+1.5, by = 2), pos=2, par("usr")[3]-0.2, xpd = TRUE, labels = rownames(toPlot), srt = 45);
            
           
        segments(mp, t(toPlot) - t(sds), mp, t(toPlot) + t(sds), lwd = 2 );
        
      } 
      else {
        barplot(toPlot, legend.text = rownames(toPlot), 
                main = paste("m/z = ", rownames(toPlot)[1], "\n", "RT = ", isoLabelReport$rt[[i]][1], "\n", "p = ", pval));
      }
    }
  }
  } # notUse
  
  ################################################
  
  dev.off();
}
