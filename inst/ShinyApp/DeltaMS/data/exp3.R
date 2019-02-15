## Experiment 3: Isotopologue guided data comparison
### Progress Bar

withProgress(message = "Analysis progress:", value  = 0, {

#Peak-picking and retention-time alignment with XCMS

  if (isoD@methodPD == "centWave") {

    incProgress(1/15, detail = "Peakdetection (centWave)")
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

    incProgress(1/15, detail = "Peakdetection (matchedFilter)")
    Sys.sleep(0.2)

    xset <- xcmsSet(files = isoD@MSFiles,
                    profmethod = isoD@profmethod,
                    polarity = isoD@polarity,
                    method = isoD@methodPD,
                    fwhm = isoD@fwhm,
                    sigma = (isoD@fwhm/2.3548),
                    max = isoD@max,
                    snthresh = isoD@snthreshMatched,
                    step = isoD@step,
                    steps = isoD@steps,
                    mzdiff = isoD@mzdiffMatched)}

### Rearrangment of xset regarding cl1, cl2, condition 1 and 2

xsetSplit<-split(xset, f=sampclass(xset))
sampclass(xsetSplit[[isoD@uC]])<-"ulSamples"
sampclass(xsetSplit[[isoD@uT]])<-"ulSamples"
sampclass(xsetSplit[[isoD@lC]])<-"lSamples"
sampclass(xsetSplit[[isoD@lT]])<-"lSamples"

xsetUpdate<-c(xsetSplit[[isoD@uC]],xsetSplit[[isoD@uT]],xsetSplit[[isoD@lC]],xsetSplit[[isoD@lT]]) #control first

### First Grouping

incProgress(1/15, detail = "First grouping")
Sys.sleep(0.2)

  xset2 <- group(xsetUpdate,
                bw=isoD@bw,
                mzwid=isoD@mzwid,
                minfrac = isoD@minfrac,
                minsamp = isoD@minsamp,
                max = isoD@maxN)


### Retention time correction
 if(isoD@plottype != "none"){ 
   pdf(file = file.path(resultsDir,"Retention time correction.pdf"))}
  if (isoD@retMeth == "obiwarp") {

    incProgress(1/15, detail = "Retention time correction (obiwarp)")
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

      incProgress(1/15, detail = "Retention time correction (Peakgroups)")
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
  dev.off()} # end of pdf retention time correction


  incProgress(1/15, detail = "Second grouping")
  Sys.sleep(0.2)

  xset2 <- group(xset2,
                 bw = isoD@bw,
                 mzwid = isoD@mzwid,
                 minfrac = isoD@minfrac,
                 minsamp = isoD@minsamp,
                 max = isoD@max)

  ### Fillpeaks

  incProgress(1/15, detail = "Fillpeaks")
  Sys.sleep(0.2)


  xset3 <- fillPeaks(xset2)

  pT<-peakTable(xset3)
  write.csv(pT, file = file.path(resultsDir, paste0(isoD@filename, ".csv")))

# Setting variables for X13CMS
  namesFolder<-c(isoD@uC,isoD@lC,isoD@uT,isoD@lT)
  folderpaths<-as.list(file.path(isoD@WD, namesFolder))
  filesInfolders<-lapply(folderpaths, function(i) {list.files(i, pattern = isoD@MSfilePattern, ignore.case = TRUE, full.names = FALSE)})
  MSpattern<-paste0(".", tail(unlist(strsplit(filesInfolders[[1]][1],"[[:punct:]]")),1))
  filesInfolders<-lapply(filesInfolders, function(i) {sub(MSpattern, "",i)})
  names(filesInfolders)<-c("unlabeled Control", "labeled Control", "unlabeled Treatment", "labeled Treatment")


sNctrl <- c(filesInfolders$`unlabeled Control`,filesInfolders$`labeled Control`)  # control samples (3 unlabeled, 3 labeled)
sNperturbation <- c(filesInfolders$`unlabeled Treatment`,filesInfolders$`labeled Treatment`) # LPS-treated samples (3 unlabeled, 3 labeled)

# labeling report for control samples:

incProgress(1/15, detail = paste("Get IsoLabelReport", isoD@condition1))
Sys.sleep(0.2)


labelsCtrl <<- getIsoLabelReport(xcmsSet = xset3,
                                 sampleNames = sNctrl,
                                 unlabeledSamples = "ulSamples",
                                 labeledSamples = "lSamples",
                                 isotopeMassDiff = (isoD@isotopeMassDiff/isoD@charge),
                                 RTwindow = isoD@Rtwindow,
                                 ppm = isoD@ppmw,
                                 massOfLabeledAtom = isoD@massOfLabeledAtom,
                                 noiseCutoff = isoD@noiseCutoff,
                                 intChoice = isoD@intChoice,
                                 varEq = isoD@varEQ,
                                 alpha = isoD@alpha,
                                 singleSample = isoD@singleSample,
                                 compareOnlyDistros = isoD@compareOnlyDistros,
                                 monotonicityTol = isoD@monoTol,
                                 enrichTol = isoD@enriTol)

# labeling report for LPS-treated samples:

incProgress(1/15, detail = paste("Get IsoLabelReport", isoD@condition2))
Sys.sleep(0.2)

labelsPerturbation <<- getIsoLabelReport(xcmsSet = xset3,
                                         sampleNames = sNperturbation,
                                         unlabeledSamples = "ulSamples",
                                         labeledSamples = "lSamples",
                                         isotopeMassDiff = (isoD@isotopeMassDiff/isoD@charge),
                                         RTwindow = isoD@Rtwindow,
                                         ppm = isoD@ppmw,
                                         massOfLabeledAtom = isoD@massOfLabeledAtom,
                                         noiseCutoff = isoD@noiseCutoff,
                                         intChoice = isoD@intChoice,
                                         varEq = isoD@varEQ,
                                         alpha = isoD@alpha,
                                         singleSample = isoD@singleSample,
                                         compareOnlyDistros = isoD@compareOnlyDistros,
                                         monotonicityTol = isoD@monoTol,
                                         enrichTol = isoD@enriTol)

# in each of the sN variables, the first 3 samples listed are of the "C12" or unlabeled type while the next 3 are of the "C13" type

classes1 <- c(rep("ulSamples", length(filesInfolders$`unlabeled Control`)), rep("lSamples",length(filesInfolders$`labeled Control`)))
classes2 <- c(rep("ulSamples", length(filesInfolders$`unlabeled Treatment`)), rep("lSamples",length(filesInfolders$`labeled Treatment`)))


# diffReport comparing labeling patterns in control vs LPS-treated samples:

incProgress(1/15, detail = "Get IsoDiffReport")
Sys.sleep(0.2)

isoDiff = getIsoDiffReport(labelsData1 = labelsCtrl,
                           labelsData2 = labelsPerturbation,
                           condition1 = isoD@condition1,
                           condition2 = isoD@condition2,
                           classes1 = classes1,
                           classes2 = classes2,
                           labeledSamples = "lSamples",
                           varEq = isoD@varEQ,
                           singleSample = isoD@singleSample)

# print labeling report to a text file (recommended to open in Excel)

incProgress(1/15, detail = "Print IsoList")
Sys.sleep(0.2)


printIsoListOutputs(listReport = labelsCtrl,
                    outputfile = file.path(resultsDir,"labelsCtrl.txt"))

# print pdf of isotopologue groups in a single labeling report plotted as relative intensity distributions
incProgress(1/15, detail = "Plot LabelReport (relative intensities)")
Sys.sleep(0.2)

plotLabelReport(isoLabelReport = labelsCtrl,
                intOption = "rel",
                classes = classes1,
                labeledSamples = "lSamples",
                outputfile = file.path(resultsDir, paste0(isoD@pdfRel,".pdf")))

# print pdf of isotopologue groups in a single labeling report plotted as absolute intensity distributions
incProgress(1/15, detail = "Plot LabelReport (absolute intensities)")
Sys.sleep(0.2)

plotLabelReport(isoLabelReport = labelsCtrl,
                intOption = "abs",
                classes = classes1,
                labeledSamples = "lSamples",
                outputfile = file.path(resultsDir, paste0(isoD@pdfAbs,".pdf")))

# print pdf of isotopologue groups from an isoDiff report comparing two conditions; plots are of relative intensity distributions
incProgress(1/15, detail = "Plot IsoDiffReport")
Sys.sleep(0.2)

plotIsoDiffReport(isoDiffReport = isoDiff,
                  xcmsSet = xset3,
                  intChoice = isoD@intChoice,
                  sampleNames1 = sNctrl,
                  sampleNames2 = sNperturbation,
                  labelReport1 = labelsCtrl,
                  labelReport2 = labelsPerturbation,
                  classes1 = classes1,
                  classes2 = classes2,
                  labeledSamples = "lSamples",
                  isotopeMassDifference = (isoD@isotopeMassDiff/isoD@charge),
                  outputfile = file.path(resultsDir,"isoDiffPlots.pdf"))

# print pdf of isotopologue groups from an isoDiff report comparing two conditions; plots are of absolute intensity distributions combined into total pools;
# reported p-values are from t-test of total pools (sum of all isotopologues in a group) between conditions
incProgress(1/15, detail = "Plot TotalIsoPools")
Sys.sleep(0.2)

plotTotalIsoPools(isoDiffReport = isoDiff,
                  xcmsSet = xset3,
                  intChoice = isoD@intChoice,
                  sampleNames1 = sNctrl,
                  sampleNames2 = sNperturbation,
                  labelReport1 = labelsCtrl,
                  labelReport2 = labelsPerturbation,
                  classes1 = classes1,
                  classes2 = classes2,
                  labeledSamples = "lSamples",
                  outputfile = file.path(resultsDir,"totalPools.pdf"))

# obtain a condensed version of the XCMS diffReport using just the C12 samples of control and LPS-treated samples; output is a matrix

incProgress(1/15, detail = "Print miniDiffReport")
Sys.sleep(0.2)

miniDiff <- miniDiffReport(xcmsSet = xset3,
                           class1sampNames = filesInfolders$`unlabeled Control`,
                           class2sampNames = filesInfolders$`unlabeled Treatment`,
                           varEq = isoD@varEQ, intChoice = isoD@intChoice)
write.csv(miniDiff, file = file.path(resultsDir,"miniDiffReport.csv"))
# filters isoDiff report to report back only isotopologue groups that are different between sample classes

incProgress(1/15, detail = "Generate filteredIsoDiff")
Sys.sleep(0.2)

filteredIsoDiff <- filterIsoDiffReport(isoDiffReport = isoDiff, alpha = isoD@alpha)

print("DeltaMS completed analysis")

})
