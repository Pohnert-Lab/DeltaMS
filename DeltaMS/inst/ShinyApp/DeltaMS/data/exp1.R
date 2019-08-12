## Experiment 1: Isotope signature
### Progress Bar

withProgress(message = "Analysis progress:", value  = 0, {

  ### Calculate list with possible isotopologues regarding iRatio

  incProgress(1/12, detail = "Get isotope ratios")
  Sys.sleep(0.2)

  AllComb <- vector("list", isoD@numAtom)
  for(i in 2:(isoD@numAtom + 1)){

    AllComb[[i]] <- IsoRatioCalc(iRatioPerc = isoD@iRatio/sum(isoD@iRatio), m = i)
  }
  AllComb[[1]] <- isoD@iRatio/sum(isoD@iRatio)

  ### Start of XCMS Procedure

  ## Peakdetection


  if (isoD@methodPD == "centWave") {

    incProgress(1/12, detail = "Peakdetection (centWave)")
    Sys.sleep(0.2)

    if(isoD@integrate == "Use filtered data"){
      integrate<-1}else{
        integrate<-2}

    xset <- xcmsSet(files = isoD@MSFiles,
                    profmethod = isoD@profmethod,
                    polarity = isoD@polarity,
                    method = isoD@methodPD,
                    ppm = isoD@ppm,
                    peakwidth=c(isoD@peakwidthMIN,isoD@peakwidthMAX),
                    snthresh = isoD@snthreshCent,
                    prefilter=c(isoD@prefilterPeakNr,isoD@prefilterIntensity),
                    mzCenterFun=isoD@mzCenterFun,
                    integrate = integrate,
                    mzdiff = isoD@mzdiffCent,
                    fitgauss = isoD@fitgauss,
                    noise = isoD@noise,
                    verbose.columns = isoD@verboseColumns)

  } else {
    incProgress(1/12, detail = "Peakdetection (matchedFilter)")
    xset <- xcmsSet(files = isoD@MSFiles,
                    profmethod = isoD@profmethod,
                    polarity = isoD@polarity,
                    method = isoD@methodPD,
                    fwhm = as.numeric(isoD@fwhm),
                    sigma = as.numeric((isoD@fwhm/2.3548)),
                    max = as.numeric(isoD@max),
                    snthresh = as.numeric(isoD@snthreshMatched),
                    step = as.numeric(isoD@step),
                    steps = as.numeric(isoD@steps),
                    mzdiff = as.numeric(isoD@mzdiffMatched))}


### In case of experiment 1 (only one folder with files)


    xset1a <- xset
    xset1b <- xset
    aSampNames <- paste0("t_", sampnames(xset))
    sampnames(xset1b)<-aSampNames
    sampclass(xset1a)<-isoD@cl1
    sampclass(xset1b)<-isoD@cl2
    xset<-c(xset1a, xset1b)



### First Grouping
    incProgress(1/12, detail = "First grouping")
  xset2 <- group(xset,
                bw=isoD@bw,
                mzwid=isoD@mzwid,
                minfrac = isoD@minfrac,
                minsamp = isoD@minsamp,
                max = isoD@maxN)


  if (isoD@NrFi!=1){
   
    if(isoD@plottype != "none"){ 
     pdf(file = file.path(resultsDir,"Retention time correction.pdf"))}
    if (isoD@retMeth == "obiwarp") {
      incProgress(1/12, detail = "Retention time correction (obiwarp)")
      if (isoD@center == 0){
        xset2 <- retcor(xset2,
                        method = isoD@retMeth,
                        plottype = isoD@plottype,
                        response = isoD@response,
                        profStep = isoD@profStep,
                        distFunc = isoD@distFunc,
                        gapInit = as.numeric(isoD@gapInit),
                        gapExtend = as.numeric(isoD@gapExtend),
                        factorDiag = as.numeric(isoD@factorDiag),
                        factorGap = as.numeric(isoD@factorGap),
                        localAlignment = as.numeric(isoD@localAlignment),
                        initPenalty = as.numeric(isoD@initPenalty))
      }else{
        xset2 <- retcor(xset2,
                        method = isoD@retMeth,
                        plottype = isoD@plottype,
                        response = isoD@response,
                        profStep = isoD@profStep,
                        distFunc = isoD@distFunc,
                        center = as.numeric(isoD@center),
                        gapInit = as.numeric(isoD@gapInit),
                        gapExtend = as.numeric(isoD@gapExtend),
                        factorDiag = as.numeric(isoD@factorDiag),
                        factorGap = as.numeric(isoD@factorGap),
                        localAlignment = as.numeric(isoD@localAlignment),
                        initPenalty = as.numeric(isoD@initPenalty))
      }


    } else {

      incProgress(1/12, detail = "Retention time correction (Peakgroups)")
      Sys.sleep(0.2)
      xset2 <- retcor(xset2,
                      method = isoD@retMeth,
                      plottype = isoD@plottypePeakgroups,
                      missing = as.numeric(isoD@missing),
                      extra = as.numeric(isoD@extra),
                      smooth = isoD@smooth,
                      span = as.numeric(isoD@span),
                      family = isoD@family)
    }
    if(isoD@plottype != "none"){ 
    dev.off()} # Finish pdf for retention time correction

    incProgress(1/12, detail = "Second grouping")
    Sys.sleep(0.2)
    xset2 <- group(xset2,
                   bw = isoD@bw,
                   mzwid = isoD@mzwid,
                   minfrac = isoD@minfrac,
                   minsamp = isoD@minsamp,
                   max = isoD@max)

  }else{
    isoD@singleSample <- TRUE
  }


  ### Fillpeaks
  incProgress(1/12, detail = "Fillpeaks")
  xset3 <- fillPeaks(xset2)
  ## Assignment of "unlabeled" and "labeled" sample name

  isoD@LabeledS   <-  isoD@cl2
  isoD@unLabeledS <- isoD@cl1

  isoD@sNctrl <-rownames(xset3@phenoData)

  ### peakTable
  xset4<-split(xset3, sampclass(xset3))
  xset5<-xset4[[1]]
  if(isoD@NrFi>1){
    xset5<-group(xset5)
  }
  sampclass(xset5)<-paste(isoD@cl1, isoD@cl2, sep = "-")
  pT<-peakTable(xset5)
  write.csv(pT, file = file.path(resultsDir, paste0(isoD@filename, ".csv")))


  ### DeltaMS Reports
  incProgress(1/12, detail = "Get IsoLabelReport")
  Sys.sleep(0.2)
  labelsCtrl <- getIsoLabelReport_doubleLabel(xcmsSet = xset3,
                                             sampleNames = isoD@sNctrl,
                                             unlabeledSamples = isoD@unLabeledS,
                                             labeledSamples = isoD@LabeledS,
                                             isotopeMassDiff = (isoD@isotopeMassDiff/isoD@charge),
                                             RTwindow = isoD@Rtwindow,
                                             ppm = isoD@ppmw,
                                             massOfLabeledAtom = isoD@massOfLabeledAtom,
                                             noiseCutoff = isoD@noiseCutoff,
                                             intChoice = isoD@intChoice,
                                             varEq = isoD@varEQ,
                                             alpha = isoD@alpha,
                                             singleSample = TRUE,
                                             compareOnlyDistros = TRUE,
                                             monotonicityTol = isoD@monoTol,
                                             enrichTol = isoD@enriTol)



  if (length(labelsCtrl$rt)>0){

    # Classify and summarize of the results

    incProgress(1/12, detail = "Get InfoList")
    Sys.sleep(0.2)
    infoList <- InfoIsoList_doubleLabel(listReport = labelsCtrl,
                                        iRatio = AllComb,
                                        errRatio = (isoD@errRatio/100),
                                        maxSD = (isoD@maxSD/100),
                                        disPeak = (isoD@isotopeMassDiff/isoD@charge),
                                        maxNumDelta = isoD@dpeak)


#     printIsoListOutputs(listReport = labelsCtrl,
#                         outputfile = file.path(resultsDir, paste0(isoD@txtReport,".txt")))
# print("IsoListOutputs done")

    incProgress(1/12, detail = "Plot LabelReport (relative intensities)")
    Sys.sleep(0.2)
    plotLabelReport_doubleLabel(isoLabelReport = labelsCtrl,
                                intOption = "rel",
                                iRatio = AllComb,
                                errRatio = isoD@errRatio,
                                disPeak = (isoD@isotopeMassDiff/isoD@charge),
                                maxNumDelta = isoD@dpeak,
                                maxLT = isoD@maxLT,
                                infoList = infoList,
                                classes = xset3@phenoData,
                                labeledSamples = isoD@LabeledS,
                                outputfile = file.path(resultsDir, paste0(isoD@pdfRel,".pdf")),
                                labIso = c(isoD@cl1,isoD@cl2))

    incProgress(1/12, detail = "Plot LabelReport (absolute intensities)")
    Sys.sleep(0.2)
    plotLabelReport_doubleLabel(isoLabelReport = labelsCtrl,
                                intOption = "abs",
                                iRatio = AllComb,
                                errRatio = (isoD@errRatio/100),
                                disPeak = (isoD@isotopeMassDiff/isoD@charge),
                                maxNumDelta = isoD@dpeak,
                                maxLT = isoD@maxLT,
                                infoList = infoList,
                                classes = xset3@phenoData,
                                labeledSamples = isoD@LabeledS,
                                outputfile = file.path(resultsDir, paste0(isoD@pdfAbs,".pdf")),
                                labIso = c(isoD@cl1,isoD@cl2))

    incProgress(1/12, detail = "Plot mass spectra")
    Sys.sleep(0.2)
    plotMassSpectrum_doubleLabel(isoLabelReport = labelsCtrl,
                                 rawFiles = isoD@MSFiles,
                                 iRatio = AllComb,
                                 errRatio = (isoD@errRatio/100),
                                 disPeak = (isoD@isotopeMassDiff/isoD@charge),
                                 maxNumDelta = isoD@dpeak,
                                 RTwin = isoD@Rtwindow ,
                                 maxLT = isoD@maxLT,
                                 infoList = infoList,
                                 classes = xset3@phenoData,
                                 labeledSamples = isoD@LabeledS,
                                 outputfile = file.path(resultsDir, paste0(isoD@pdfMZ,".pdf")),
                                 labIso = c(isoD@cl1,isoD@cl2))


    incProgress(1/12, detail = "Plot extracted ion current chromatograms")
    Sys.sleep(0.2)
    plotEIC_doubleLabel(isoLabelReport = labelsCtrl,
                        rawFiles = isoD@MSFiles,
                        iRatio = AllComb,
                        errRatio = (isoD@errRatio/100),
                        disPeak = (isoD@isotopeMassDiff/isoD@charge),
                        maxNumDelta = isoD@dpeak,
                        RTwin = isoD@Rtwindow,
                        maxLT = isoD@maxLT,
                        infoList = infoList,
                        classes = xset3@phenoData,
                        labeledSamples = isoD@LabeledS,
                        outputfile =  file.path(resultsDir, paste0(isoD@pdfEIC,".pdf")),
                        labIso = c(isoD@cl1,isoD@cl2))


  }

#Last Function that is called: Write "Analysis settings.html

 print("DeltaMS completed analysis")
  })
