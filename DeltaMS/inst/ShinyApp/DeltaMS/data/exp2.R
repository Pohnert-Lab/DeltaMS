## Experiment 2: Dataset comparison

### Progress Bar
withProgress(message = "Analysis progress:", value  = 0, {

#Peak-picking and retention-time alignment with XCMS

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
  Sys.sleep(0.2)

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


### First Grouping

incProgress(1/12, detail = "First grouping")
Sys.sleep(0.2)

xset2 <- group(xset,
                bw=isoD@bw,
                mzwid=isoD@mzwid,
                minfrac = isoD@minfrac,
                minsamp = isoD@minsamp,
                max = isoD@maxN)


### Retention time correction
if(isoD@plottype != "none"){ 
pdf(file = file.path(resultsDir,"Retention time correction.pdf"))}

if (isoD@retMeth == "obiwarp") {
  incProgress(1/12, detail = "Retention time correction (obiwarp)")
  Sys.sleep(0.2)

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
  }} else {
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
dev.off()} #end of pdf retention time correction

incProgress(1/12, detail = "Second grouping")
Sys.sleep(0.2)

xset2 <- group(xset2,
                bw = isoD@bw,
                mzwid = isoD@mzwid,
                minfrac = isoD@minfrac,
                minsamp = isoD@minsamp,
                max = isoD@max)

### Fillpeaks
incProgress(1/12, detail = "Fillpeaks")
Sys.sleep(0.2)


xset3 <- fillPeaks(xset2)

## Assignment of "unlabeled" and "labeled" sample name
pT<-peakTable(xset3)
write.csv(pT, file = file.path(resultsDir, paste0(isoD@filename, ".csv")))

isoD@sNctrl <-rownames(xset3@phenoData)

### DeltaMS Reports
incProgress(1/12, detail = "Get LabelReport")
Sys.sleep(0.2)

labelsCtrl <- getIsoLabelReportdeltaMS(
  xcmsSet = xset3,
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
  singleSample = FALSE,
  compareOnlyDistros = TRUE,
  monotonicityTol = isoD@monoTol,
  enrichTol = isoD@enriTol
)

if (length(labelsCtrl$rt)>0){

  # Classify and summarize of the results
  incProgress(1/12, detail = "Get InfoList")
  Sys.sleep(0.2)


  infoList <- InfoIsoList(listReport = labelsCtrl,
                                      iRatio = isoD@iRatio,
                                      errRatio = (isoD@errRatio/100),
                                      maxSD = (isoD@maxSD/100),
                                      disPeak = (isoD@isotopeMassDiff/isoD@charge),
                                      maxNumDelta= isoD@dpeak)

  # printIsoListOutputs(listReport = labelsCtrl,
  #                     outputfile = file.path(resultsDir, isoD@txtReport))


  incProgress(1/12, detail = "Plot LabelReport (relative intensities)")
  Sys.sleep(0.2)

  plotLabelReport_isoRat(isoLabelReport = labelsCtrl,
                              intOption = "rel",
                              iRatio = isoD@iRatio,
                              errRatio = (isoD@errRatio/100),
                              disPeak = (isoD@isotopeMassDiff/isoD@charge),
                              maxNumDelta = isoD@dpeak,
                              maxLT = isoD@maxLT,
                              infoList = infoList,
                              classes = xset3@phenoData,
                              labeledSamples = isoD@LabeledS,
                              outputfile = file.path(resultsDir, paste0(isoD@pdfRel,".pdf")),
                              labIso = c(isoD@cl2,isoD@cl1))

  incProgress(1/12, detail = "Plot LabelReport (absolute intensities)")
  Sys.sleep(0.2)

  plotLabelReport_isoRat(isoLabelReport = labelsCtrl,
                              intOption = "abs",
                              iRatio = isoD@iRatio,
                              errRatio = (isoD@errRatio/100),
                              disPeak = (isoD@isotopeMassDiff/isoD@charge),
                              maxNumDelta = isoD@dpeak,
                              maxLT =isoD@maxLT,
                              infoList = infoList,
                              classes = xset3@phenoData,
                              labeledSamples = isoD@LabeledS,
                              outputfile = file.path(resultsDir, paste0(isoD@pdfAbs,".pdf")),
                              labIso = c(isoD@cl2,isoD@cl1))


  incProgress(1/12, detail = "Plot mass spectra")
  Sys.sleep(0.2)

  plotMassSpectrum(isoLabelReport = labelsCtrl,
                               rawFiles = isoD@MSFiles,
                               iRatio = isoD@iRatio,
                               errRatio = (isoD@errRatio/100),
                               disPeak = (isoD@isotopeMassDiff/isoD@charge),
                               maxNumDelta = isoD@dpeak,
                               RTwin = isoD@Rtwindow ,
                               maxLT = isoD@maxLT,
                               infoList = infoList,
                               classes = xset3@phenoData,
                               labeledSamples = isoD@LabeledS,
                               outputfile = file.path(resultsDir, paste0(isoD@pdfMZ,".pdf")),
                               labIso = c(isoD@cl2,isoD@cl1))


  #plotLabelReport(isoLabelReport = labelsCtrl, intOption = "rel", classes = xset3@phenoData, labeledSamples = LabeledS, outputfile = pdfRel)

  incProgress(1/12, detail = "Plot extracted ion current chromatograms")
  Sys.sleep(0.2)

  plotExtractedIonChromatogram(isoLabelReport = labelsCtrl,
                      rawFiles = isoD@MSFiles,
                      iRatio = isoD@iRatio,
                      errRatio = (isoD@errRatio/100),
                      disPeak = (isoD@isotopeMassDiff/isoD@charge),
                      maxNumDelta = isoD@dpeak,
                      RTwin = isoD@Rtwindow,
                      maxLT = isoD@maxLT,
                      infoList = infoList,
                      classes = xset3@phenoData,
                      labeledSamples = isoD@LabeledS,
                      outputfile =  file.path(resultsDir, paste0(isoD@pdfEIC,".pdf")),
                      labIso = c(isoD@cl2,isoD@cl1))
}

print("DeltaMS completed analysis")

})
