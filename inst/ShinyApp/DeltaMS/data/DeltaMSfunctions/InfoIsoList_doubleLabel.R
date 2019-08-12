
# Summary of getIsoLabelReport() of a double labeled sample
InfoIsoList_doubleLabel <- function(listReport, iRatio, errRatio, disPeak, maxNumDelta, maxSD) {
  # isotopic ratio # 18.01.2016
  if (is.null(iRatio))  {
    iRatio <- c(0.5, 0.5);
  } 
  if (is.null(errRatio))  {
    errRatio <-0.05;
  } 
  # Permissible standard deviation within a group. # August 2017
  if (is.null(maxSD))  {
    maxSD <-0.1;
  } 
  if (maxNumDelta < 1) {
    maxNumDelta <- 1
  }
  
  # 
  numGroups = length(listReport[[1]]); # Number of ISO - groups
  # Number of isotopic pattern 
  numPat <- length(iRatio)
  
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
  
  curInt <- min(c(interMax, maxNumDelta))

# Initialization of the output variable
  iMat <- matrix(data = NA, nrow = numGroups, ncol = 5)
  cNames <- c("distance12", "numPeaks", "found", "maxdiff", "rank")
  colnames(iMat)<-cNames
   
   
  iMatDist<- matrix(data = NA, nrow = numGroups, ncol = curInt)
  iMatPat <- matrix(data = NA, nrow = numGroups, ncol = (curInt-0)) # -1
  
  curPat <- curInt #- 1
  
  # 	interval boundaries
  for (ic in 1:numPat) {
    if (ic > curPat) {
      break
    }
    B <- iRatio[[ic]]
    BlI <- B*(1-errRatio)
    BuI <- B*(1+errRatio)
    numTest <- length(B)
     
    bt <-B

    for (i in 1:numGroups) {
      # Ratio between . peaks
      toPlot = cbind(listReport$meanRelU[[i]]);
      ratA <- toPlot/sum(toPlot)
      numPeaks = dim(toPlot)[1]
#     mean value      
#      sdmean <- mean(cbind(listReport$sdRelU[[i]]))  # groups with high SD should be discarded
# max value      
      sdmean <- max(cbind(listReport$sdRelU[[i]]))
      if(is.na(sdmean)){
        sdmean <- 0
      }
 #     if (i == 54) {
#        print(i)
#      }
      if (numPeaks < numTest){
         next
      }
      diffA <- diff(listReport$isotopologue[[i]])
      isDiff <- round(diffA[1])
   
      idxPeak <- (round((listReport$isotopologue[[i]]-listReport$isotopologue[[i]][1])/disPeak))+1
      isDiffA <- round(diffA/disPeak)

      expRatio <- 1
       curTest <- min(c(numPeaks, numTest))
      for (ip in 1:curTest){
         if (!(ratA[ip] >= BlI[ip] && ratA[ip] <= BuI[ip])) {
            expRatio <- 0
         }
      }
      # test: if correct mz-distance
      if (expRatio == 1){
        for (ip in 1:(curTest-1)){
           if ((isDiffA[ip] > 1)){
             expRatio <- 0
           }
        }
      }
 
      iMatPat[i, ic]<-expRatio
      if (ic == 1) {
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
      if ((length(idxPeak) > 1) && (max(diff(idxPeak))==1)) {
        iMat[i, 4] <- 1
      } else {
        iMat[i, 4] <- 0
      }
      if (iMat[i,1]==1  && iMat[i,2]>2 && iMat[i,4] == 0) {
        zz <- 0
        for (k in 1:min(iMat[i,2], curInt)) {  # Error: Subscript out or range (27 june 2017)
         # print(c(i, k))
          if (!(is.na(iMatDist[i, k]))){
            zz <- zz + 1
          } else {
            break
          }
        }  
       iMat[i, 5] <- zz
      } else {
        iMat[i, 5] <- 0
      }
    } 
  }
  infoList = list( iDist = iMatDist, iPat = iMatPat, infoMat = iMat)

  return(infoList);  
}

