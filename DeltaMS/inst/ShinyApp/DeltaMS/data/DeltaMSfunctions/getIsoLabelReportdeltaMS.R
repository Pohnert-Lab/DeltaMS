# outputs report of features enriched for labeled isotope in samples named in "sampleNames"
getIsoLabelReportdeltaMS <- function(xcmsSet, sampleNames, unlabeledSamples, labeledSamples, isotopeMassDiff, RTwindow, ppm, 
                              massOfLabeledAtom, noiseCutoff, intChoice = "intb", varEq = FALSE, alpha, singleSample = FALSE, compareOnlyDistros = FALSE,
                              monotonicityTol = FALSE, enrichTol = 0.1) {
  # Changes from original 
  #  Error handling, if empty lists are the result.
  #
  # Variables:
  # xcmsSet: xcmsSet object containing grouped and retention-time-aligned peaks (i.e. after calling group() and retcor() in XCMS); 
  #          sample classes should be designated as "unlabeledSamples" and "labeledSamples";
  #          may contain samples from one or more biological conditions
  # sampleNames: vector of names of the unlabeled and labeled samples in the condition for which the report is being generated;
  #              taken from rownames of xcmsSet@phenoData
  # unlabeledSamples: character variable name of unlabeled samples as designated in xcmsSet@phenoData (e.g. "C12")
  # labeledSamples: character variable name of labeled samples (e.g. "C13")
  # isotopeMassDiff: difference in mass between labeled and unlabeled atom (e.g. 1.00335 for C13) 
  # RTwindow: retention time window in which all peaks are considered to be co-eluting
  # ppm: ppm allowance for deviation of peaks within an isotopologue group from expected m/z; instrument dependent
  # massOfLabeledAtom: e.g. 12.0000 for C12
  # noiseCutoff: ion intensity cutoff below which a peak is considered noise; instrument dependent
  # intChoice: choice of "maxo", "into", or "intb" for peak intensities
  # varEq: Boolean indicating whether to assume equal variances for the peak intensities used in the t-test comparing unlabeled
  #        vs labeled samples
  # alpha: p-value cutoff for significance of label enrichment
  # singleSample: Boolean indicating whether only a single sample exists for the unlabeled and labeled conditions; 
  #               if TRUE, all potential isotopologue groups are returned and no statistics are performed
  # compareOnlyDistros: Boolean indicating whether to assess for true labeling patterns (FALSE) or 
  #                     to compare label distributions if both samples are labeled (TRUE)
  # monotonicityTol: tolerance parameter used to enforce expected ion intensity pattern (i.e. monotonic decrease from M0 to Mn) 
  #                  in unlabeled samples; a low value closer to 0 enforces stricter monotonicity; default is to not enforce monotonicity
  #                  due to potential carryover between samples
  # enrichTol: tolerance parameter for enforcing enrichment of higher isotopologues in labeled samples; a value of 0 enforces strict requirement
  #            for enrichment of higher isotopologues to be higher in labeled samples
  
  ### UNPACK THE XCMS OBJECT ###
  groups = data.frame(xcmsSet@groups);
  peakIntensities = getPeakIntensities(xcmsSet, sampleNames, intChoice);
  peakIntensities = peakIntensities[order(groups$mzmed),];
  groups = groups[order(groups$mzmed),]; # order peaks by m/z
  groupRTs = groups$rtmed;
  groupMzs = groups$mzmed;
  groupIDs = as.numeric(rownames(groups));
  nGroups = length(groupMzs);
  classes = as.character(xcmsSet@phenoData[match(sampleNames,rownames(xcmsSet@phenoData)),]);
  numSamples = length(classes);
  
  intensities1= peakIntensities[,which(classes == unlabeledSamples), drop = FALSE];
  intensities2 = peakIntensities[,which(classes == labeledSamples), drop = FALSE];
  
  iMD = isotopeMassDiff; 
  
  ### ISOTOPOLOGUE SEARCH: ###
  # temporary storage lists for potential isotopologue pairs
  base = list(); # index of base peak
  labeled = list(); # index of labeled peak
  basePeak = list(); # m/z of base peak
  labeledPeak = list(); # m/z of labeled peak
  
  groupIndicesByRT = order(groupRTs);
  orderedGroupRTs = groupRTs[groupIndicesByRT];
  # the search:
  for (i in 1:nGroups) {
    binI = groupIndicesByRT[orderedGroupRTs - orderedGroupRTs[i] >= 0 & orderedGroupRTs - orderedGroupRTs[i] <= RTwindow];
    bin = groups[binI,];
    binSize = length(binI);
    I = groupIndicesByRT[i];
    if (binSize > 0) {
      # do pairwise comparisons between every m/z group represented in that bin
      # a = unlabeled peak
      # b = labeled peak
      for (j in 1:binSize) {
        if (groups$mzmed[I] < bin$mzmed[j]) {
          a = I;
          b = binI[j];
        }
        else {
          a = binI[j];
          b = I;
        }   
        
        
   
        delta = (groupMzs[b] - groupMzs[a])/iMD;
        DELTA = round(delta); # candidate number of labeled atoms
        if (DELTA == 0) {
          next;
        }
        # if delta is a multiple of mass defect to within ppm error:
        if ( delta <= DELTA*(1+ppm/1000000) + (groupMzs[a]*ppm/1000000)/(iMD*(1-ppm/1000000)) 
             && delta >= DELTA*(1-ppm/1000000) - (groupMzs[a]*ppm/1000000)/(iMD*(1+ppm/1000000)) ) {
          # ...check if mass defect is too large
          if ( DELTA  * massOfLabeledAtom >= groupMzs[a]) {
            next;
          }
          # ...check that intensities of labeled peak are less than those of base peak in the unlabeled samples (class == 1)
          if ( mean(intensities1[b,]) > mean(intensities1[a,]) && !compareOnlyDistros ) {
            next;
          } 
          # ...check that intensities are not 0 in all samples
          if ( all(intensities1[a,] == 0) && all(intensities2[a,] == 0)) {
            next;
          }
          if ( all(intensities1[b,] == 0) && all(intensities2[b,] == 0)) {
            next;
          }
          # record pair of peaks if all filters are passed
          base = c(base, a);
          labeled = c(labeled, b);
          basePeak = c(basePeak, groupMzs[a]);
          labeledPeak = c(labeledPeak, groupMzs[b]);
        }
      }
    }    
  }
  if (length(base) == 0) {
    return (NULL)
  }
  
  labelsMatrix = as.matrix(cbind(unlist(base), 
                                 unlist(labeled), 
                                 unlist(basePeak), 
                                 unlist(labeledPeak)));
  labelsMatrix = labelsMatrix[order(labelsMatrix[,3],labelsMatrix[,4]),]; # order into isotopologue groups

  ### DATA CLEAN-UP: ###
  # Part I:
  # Remove duplicate pairs and pairs that are both labeled peaks for other base peaks
  numPutativeLabels = dim(labelsMatrix)[1];
  if (is.null(numPutativeLabels)) {
    basePeaks = labelsMatrix[1];
  } else {
    basePeaks = unique(labelsMatrix[,1]);
  }
  numLabeledPeaks = length(basePeaks);
  outtakes = list();
  for ( i in 1:numPutativeLabels ) { 
    B = labelsMatrix[,2] == labelsMatrix[i,1];
    A = labelsMatrix[B,1];
    C = which(labelsMatrix[,1] %in% A);
    if ( any(labelsMatrix[C,2] == labelsMatrix[i,2]) ) {
      outtakes = c(outtakes,i);
      next;
    }
    if ( i < numPutativeLabels ) {
      A = (i+1):numPutativeLabels;
      if ( any(labelsMatrix[A,1] == labelsMatrix[i,1] && labelsMatrix[A,2] == labelsMatrix[i,2])) {
        outtakes = c(outtakes,i);
      }
    }
  }
  outtakes = unlist(outtakes);
  if (!is.null(outtakes)){
    labelsMatrix = labelsMatrix[-outtakes,];
  }

  # Part II:
  # Check for significant difference between labeling patterns in unlabeled and labeled samples;
  # record only those groups that are different at p-value = alpha cutoff
  numPutativeLabels = dim(labelsMatrix)[1];
  basePeaks = unique(labelsMatrix[,1]);
  numLabeledPeaks = length(basePeaks);

  # output lists:
  base = list(); # base peak associated with isotopologue group
  mz = list(); # m/z 
  ID = list(); # peak ID number in xcmsSet@groups; if user wants to plot EICs of isotopologues, input this to XCMS's getEIC routine
  RT = list(); # retention time
  absInt1 = list(); # mean absolute intensity of peak in unlabeled samples
  absInt2 = list(); # ibid in labeled samples
  relInt1 = list(); # mean relative intensity of peak in unlabeled samples
  relInt2 = list(); # ibid in labeled samples
  totInt1 = list(); # total absolute intensity of all peaks in isotopologue group in unlabeled samples
  totInt2 = list(); #ibid in labeled samples
  CVabsInt1 = list(); # coefficient of variation of total ion intensity of isotopologue group in unlabeled samples
  CVabsInt2 = list(); # ibid in labeled samples
  SDrelInt1 = list(); # std dev of relative intensity of each isotopologue peak in unlabeled samples
  SDrelInt2 = list(); # ibid in labeled samples  
  foldEnrichment = list(); # (relative intensity in labeled samples) / (relative intensity in unlabeled samples)
  pvalues = list(); # from t-test for difference between unlabeled and labeled samples 
  sampleIntensities = list(); # ion counts (of type "intChoice") for each peak in isotopologue group in every sample
  
  # process each isotopologue group:
  j = 1;
  
  for (i in 1:numLabeledPeaks) {
    a = basePeaks[i]; # groupID of base peak
    baseIntensities = c(intensities1[a,], intensities2[a,]);
    isotopologues = list();
    IDs = list();
    RTs = list();
    numisotopologues = 0;
    k = j;
    # count the number of labeled peaks in each group
    while (k <= numPutativeLabels && labelsMatrix[k,1] == a) {
      isotopologues = c(isotopologues, groupMzs[labelsMatrix[k,2]]);
      IDs = c(IDs, groupIDs[labelsMatrix[k,2]]);
      RTs = c(RTs, groupRTs[labelsMatrix[k,2]]);
      numisotopologues = numisotopologues + 1;
      k = k+1;
    }
    if (numisotopologues == 0)
    {
      next();
    }
    
    isotopologues = unlist(isotopologues);
    IDs= unlist(IDs);
    RTs = unlist(RTs); 

    # discard peaks of low intensities in unlabeled samples as candidate base peaks
    if ( mean(intensities1[a,]) < noiseCutoff) { 
      j = k;
      next;
    }
    
    # create and fill a matrix to contain intensities of labeled peaks 
    abs1 = list();
    abs2 = list();
    
    labeledIntensities = matrix(rep(0,numisotopologues*numSamples),nrow = numisotopologues, ncol = numSamples);
    
    for (l in 1:numisotopologues) {
      b = labelsMatrix[j+l-1,2]; # groupID of labeled peak
      labeledIntensities[l,] = cbind(intensities1[b,,drop=FALSE], intensities2[b,,drop=FALSE]);
      abs1 = c(abs1, mean(intensities1[b,]));
      abs2 = c(abs2, mean(intensities2[b,]));
    }
    abs1 = unlist(abs1);
    abs2 = unlist(abs2);

    # remove redundantly called isotopologues by keeping only those with the smallest ppm error from expected m/z
    if (numisotopologues != length(unique(round(isotopologues)))) {
      M0 = round(groupMzs[a]);
      isos = round(isotopologues);
      reduced = unique(isos);
      numUniqIsos = length(reduced);
      outtakes = list();
      for (r in 1:numUniqIsos) {
        q = which(isos == reduced[r]);
        if (length(q) > 1 ) {
          massdefect = iMD * (reduced[r] - M0);
          delta = abs(groupMzs[IDs[q]] - groupMzs[a] - massdefect);
          outtakes = c(outtakes, q[which(delta != min(delta))]);          
        }
      }
      if (length(outtakes) >0) {
        outtakes = unlist(outtakes);
        isotopologues = isotopologues[-outtakes];
        IDs = IDs[-outtakes];
        RTs = RTs[-outtakes];
        numisotopologues = length(isotopologues);
        abs1 = abs1[-outtakes];
        abs2 = abs2[-outtakes];
        labeledIntensities = labeledIntensities[-outtakes,,drop=FALSE];
      }
    }

    # remove isotopologues whose intensities do not decrease monotonically in unlabeled samples from M0 to Mx   
    if (!compareOnlyDistros && monotonicityTol) {
      meanMprevUL = mean(intensities1[a,]); # mean intensity of M0 in unlabeled samples
      outtakes = list();
      for (l in 1:numisotopologues) {
        if (l ==1 && abs1[l] > (1+monotonicityTol)*meanMprevUL) {
          outtakes = c(outtakes,l);
        }
        else if (l > 1 && abs1[l] > (1+monotonicityTol)*meanMprevUL && round(isotopologues[l] - isotopologues[l-1]) > 1) {
          outtakes = c(outtakes,l);
        }
        else {
          meanMprevUL = abs1[l];
        }
      }
      outtakes = unlist(outtakes);
      if (length(outtakes) > 0) {
        abs1 = abs1[-outtakes];
        abs2 = abs2[-outtakes];
        labeledIntensities = labeledIntensities[-outtakes,,drop=FALSE];
        isotopologues = isotopologues[-outtakes];
        IDs= IDs[-outtakes];
        RTs = RTs[-outtakes];
      }
    }
    
    isotopologues = c(groupMzs[a], isotopologues);
    IDs= c(groupIDs[a], IDs);
    RTs = c(groupRTs[a], RTs); 
    allIntensities = rbind(baseIntensities, labeledIntensities);
    abs1 = c(mean(intensities1[a,]), abs1);
    abs2 = c(mean(intensities2[a,]), abs2);

    numisotopologues = length(isotopologues);
    sumIntensities = colSums(allIntensities);
    tot1 = mean(sumIntensities[1:dim(intensities1)[2]]); 
    tot2 = mean(sumIntensities[(dim(intensities1)[2]+1):numSamples]);
    cv1 = sd(sumIntensities[1:dim(intensities1)[2]])/tot1;
    cv2 = sd(sumIntensities[(dim(intensities1)[2]+1):numSamples])/tot2;

 
    # obtain relative intensities
    groupIntensities = allIntensities/ 
      matrix(rep(sumIntensities,numisotopologues),nrow = numisotopologues, byrow = TRUE);
    
    gI1 = groupIntensities[,1:dim(intensities1)[2], drop = FALSE];
    gI2 = groupIntensities[,(dim(intensities1)[2]+1):numSamples, drop = FALSE];
    
    # discard any samples with no signal for the isotopologue group    
    gI1 = gI1[,colSums(is.na(gI1))==0, drop = FALSE];
    gI2 = gI2[,colSums(is.na(gI2))==0, drop = FALSE];

    # discard whole group if majority of samples lack signal for the group
    if (dim(gI1)[2] < dim(intensities1)[2]/2 || dim(gI2)[2] < dim(intensities2)[2]/2) {
      j = k;
      next;
    }
    
    # calculate mean relative intensities
    rel1 = rowMeans(gI1);
    rel2 = rowMeans(gI2);
    sd1 = apply(gI1, 1, sd);
    sd2 = apply(gI2, 1, sd);

    enrichRatios = rel2/rel1;
    
    # discard whole group if enrichment is higher in unlabeled samples than in labeled ones, by a factor of 1+enrichTol
    if (!compareOnlyDistros) {
      if (enrichRatios[1] > (1+enrichTol)) {
        j = k;
        next;
      }
    }
    
    if (!singleSample) {
      pvalue = list();
      for (l in 1:numisotopologues) {
        if( all(gI1[l,] == 1) && all(gI2[l,] == 0) || is.infinite(enrichRatios[l]) ) {
          pvalue = c(pvalue,0);
        }
        else {
          T = try(t.test(gI1[l,], gI2[l,], var.equal = varEq), silent = TRUE);
          if (class(T) == "try-error") {
            pvalue = c(pvalue, 1);
            break;
          }
          else {
            pvalue = c(pvalue, T$p.value);
          }
        }
      }     

      # store data if at least one labeled isotopologue is enriched in labeled samples
      if ( any(unlist(pvalue) < alpha) && !any(unlist(pvalue) == 1) ) {
        base = c(base, groupMzs[a]);
        mz = c(mz, list(isotopologues));
        ID = c(ID, list(IDs)); 
        RT = c(RT, list(RTs)); 
        absInt1 = c(absInt1, list(abs1));
        absInt2 = c(absInt2, list(abs2)); 
        relInt1 = c(relInt1, list(rel1)); 
        relInt2 = c(relInt2, list(rel2));
        CVabsInt1 = c(CVabsInt1, cv1);
        CVabsInt2 = c(CVabsInt2, cv2); 
        totInt1 = c(totInt1, tot1);
        totInt2 = c(totInt2, tot2);
        SDrelInt1 = c(SDrelInt1, list(sd1));
        SDrelInt2 = c(SDrelInt2, list(sd2)); 
        foldEnrichment = c(foldEnrichment, list(enrichRatios)); 
        pvalues = c(pvalues, list(unlist(pvalue)));
        sampleIntensities = c(sampleIntensities,list(allIntensities));
      }
    }
    else { # if only single replicate is available for unlabeled and labeled samples, record all isotopologue groups
      deltaSpec = sum(abs(rel1[1:numisotopologues-1] - rel2[1:numisotopologues-1])); # differences in relative intensity profiles
      base = c(base, groupMzs[a]);
      mz = c(mz, list(isotopologues));
      ID = c(ID, list(IDs)); 
      RT = c(RT, list(RTs)); 
      absInt1 = c(absInt1, list(abs1));
      absInt2 = c(absInt2, list(abs2)); 
      relInt1 = c(relInt1, list(rel1)); 
      relInt2 = c(relInt2, list(rel2));
      CVabsInt1 = c(CVabsInt1, cv1);
      CVabsInt2 = c(CVabsInt2, cv2);
      totInt1 = c(totInt1, tot1);
      totInt2 = c(totInt2, tot2);
      SDrelInt1 = c(SDrelInt1, list(sd1));
      SDrelInt2 = c(SDrelInt2, list(sd2));
      foldEnrichment = c(foldEnrichment, list(enrichRatios)); 
      pvalues = c(pvalues, deltaSpec); # not real p-values
      sampleIntensities = c(sampleIntensities,list(allIntensities));
    }
    j = k;
  }

  labelsData = list(compound = base, isotopologue = mz, groupID = ID, rt = RT, meanAbsU = absInt1, totalAbsU = totInt1, cvTotalU = CVabsInt1,
                    meanAbsL = absInt2, totalAbsL = totInt2, cvTotalL = CVabsInt2, meanRelU = relInt1, meanRelL = relInt2, p_value = pvalues, enrichmentLvsU = foldEnrichment,
                    sdRelU = SDrelInt1, sdRelL = SDrelInt2, sampleData = sampleIntensities);
  
  return(labelsData);  
}


