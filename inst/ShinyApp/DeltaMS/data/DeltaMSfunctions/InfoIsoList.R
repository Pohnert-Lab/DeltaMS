# Summary of getIsoLabelReport() 
InfoIsoList <- function(listReport, iRatio, errRatio, disPeak, maxNumDelta, maxSD) {
  # isotopic ratio # 18.01.2016
  if (is.null(iRatio))  {
    iRatio <- c(0.5, 0.5);
  } 
  # Permissible standard deviation within a group. # August 2017
  if (is.null(maxSD))  {
    maxSD <-0.1;
  } 
  
  if (is.null(errRatio))  {
    errRatio <-0.05;
  } 
  #
  if (maxNumDelta < 1) {
    maxNumDelta <- 1
  }
    
  
  iRatio <- iRatio/sum(iRatio)
  # 
  numGroups = length(listReport[[1]]); # Number of ISO - groups
  # Number of isotopic pattern 
  numPat <- 1 # always one if samples labled with one isotope # length(iRatio)
  
  interMax <- 1
  # Intervall of interest for the isotopologue pattern 
  for (i in 1:numGroups) {
    df <- listReport$isotopologue[[i]]
    ldf <- length(df)
    fDist <- round((df[ldf] - df[1])/disPeak)
    if (fDist > interMax){
      interMax <- fDist
    }
  }
  
  curInt <- min(c(interMax, (maxNumDelta+2)) ) # +2 Error out of bound: 2798
  
  # Initialization of the output variable
  iMat <- matrix(data = NA, nrow = numGroups, ncol = 4)
  cNames <- c("distance12", "numPeaks", "found", "IDXdiff")
  colnames(iMat)<-cNames
  
  
  iMatDist<- matrix(data = NA, nrow = numGroups, ncol = curInt)
  iMatPat <- matrix(data = NA, nrow = numGroups, ncol = (curInt-1))

  
  curPat <- curInt - 1
  
  # 	interval boundaries
  for (ic in 1:numPat) {
    if (ic > curPat) {
      break
    }
    B <- iRatio
    BlI <- B*(1-errRatio)
    BuI <- B*(1+errRatio)
    numTest <- length(B)
    
    ratA <-B
    
    for (i in 1:numGroups) {
      # Ratio between . peaks
      toPlot = cbind(listReport$meanRelU[[i]], listReport$meanRelL[[i]] );
      numPeaks = dim(toPlot)[1]
      if (numPeaks < numTest){
        next
      }
      sdmeanU <- (cbind(listReport$sdRelU[[i]]))  # groups with high SD should be discarded
      sdmeanL <- (cbind(listReport$sdRelL[[i]]))  
      sdmean <- 0
      for (j in 1:(numPeaks-1)) {
# mean value        
#       sdmean <- sdmeanU[j,1]+sdmeanL[j+1,1]
# max value        
        sdmean <- max(c(sdmean, sdmeanU[j,1], sdmeanL[j+1,1]))
      }
#      sdmean <- sdmean/(2*(numPeaks-1)) # mean value
      if(is.na(sdmean)){
        sdmean <- 0
      }
      
#      if (i == 27) {
#      print(i)
#      }
      diffA <- diff(listReport$isotopologue[[i]])
      isDiff <- round(diffA[1])
      
      idxPeak <- (round((listReport$isotopologue[[i]]-listReport$isotopologue[[i]][1])/disPeak))+1
      isDiffA <- round(diffA/disPeak)
# Error out of bounds
      numPeaks <- min(c((numPeaks), (maxNumDelta+1)))
      for (j in 1:(numPeaks-1)) {
        if ((toPlot[j, 1] + toPlot[j+1, 2]) > 0) {
           ratA[1] <- toPlot[j, 1]/(toPlot[j, 1] + toPlot[j+1, 2])
           ratA[2] <- toPlot[j+1, 2]/(toPlot[j, 1] + toPlot[j+1, 2])
        } else {
          ratA[1] <- 0
          ratA[2] <- 0
        }
        #ratA <- toPlot/sum(toPlot)
      
         expRatio <- 1
         curTest <- min(c(numPeaks, numTest))
         for (ip in 1:curTest){
           if (!(ratA[ip] >= BlI[ip] && ratA[ip] <= BuI[ip])) {
              expRatio <- 0
           }
         }
         # test: if correct mz-distance
  #     if (expRatio == 1){
  #         for (ip in 1:(curTest-1)){
  #            if ((isDiffA[j] > 1)){
  #              expRatio <- 0
  #            }
  #         }
  #       }
      # ic mit j ersetzen ?
      iMatPat[i, j]<-expRatio  # Error: Out of bounds
      if (j == 1) {
        id <- which((idxPeak) > curInt)
        if ((length(id) > 0) && (id > 0)) {
          idxPeak <- idxPeak[-id]
        }
        #iMatDist[i, 1] <- 1
        iMatDist[i, idxPeak] <- 1
        iMat[i, 1] <- round(diffA[1]/disPeak) #isDiff / 2
        iMat[i, 2] <- numPeaks
      }
      if (expRatio > 0) {
        iMat[i, 3] <- ic
      } else {
        if (is.na(iMat[i,3])) {
          iMat[i, 3] <- 0
        } 
      }
      if (sdmean > maxSD) {  # Is the permissible standard deviation within a group exceeded?
        iMat[i, 3] <- 0
      }
      if (expRatio > 0) {
        iMat[i, 4] <- j
      } else {
        if (is.na(iMat[i,4])) {
          iMat[i, 4] <- 0
        } 
      }
      if (iMat[i, 1] > maxNumDelta) {  # max. distance of peaks
        iMat[i, 3] <- 0
      }
     
#      if (iMat[i,1]==1  && iMat[i,2]>2 && iMat[i,4] == 0) {
#        zz <- 0
#        for (k in 1:iMat[i,2]) {
#          if (!(is.na(iMatDist[i, k]))){
#            zz <- zz + 1
#          } else {
#            break
#          }
#        }  
#        iMat[i, 5] <- zz
#      } else {
#        iMat[i, 5] <- 0
#      }
      } 
    }
  }
  infoList = list( iDist = iMatDist, iPat = iMatPat, infoMat = iMat)
  
  return(infoList);  
}
