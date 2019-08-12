
###############################################################
### Source file for DeltaMS containing additional functions ###
###############################################################

#######################################
#S4 class with DeltaMS variables ###
#######################################

InputDefaults<- setClass(
  "inputDefaultsDeltaMS",

  slots = c(
	#Additional variables

	NrFo				= "numeric",
	MSFiles     		= "character",
	MSfilePattern 		= "character",
	NrFi         		= "numeric",
	restore    			= "logical",
	folderNames 		= "character",
	saveFileName 		= "character",
	CorrectFilePath 	= "character",
	CorrectFileName 	= "character",
	WD           		= "character",
	foldersOld   		= "character",


  #General settings
	higherAbundIso		= "numeric",
	lowerAbundIso		= "numeric",
	profmethod  		= "character",
	profparam 			= "character",
	mslevel   			= "numeric",
	nSlaves   			= "numeric",
	TypeofExp       	= "numeric",
	methodPD   			= "character",
	isotopeMassDiff  	= "numeric",
	retMeth  			= "character",
	hp            		= "numeric",
	wp           		= "numeric",
	iOdm           		= "character",
	manOrCalc     		= "character",
	charge				= "numeric",
	isotope1      		= "character",
	isotope2      		= "character",



    ##centWave arguments
    mzdiffCent   		= "numeric",
    lockMassFreq  		= "logical",
    polarity  			= "character",
    peakwidthMIN   		= "numeric",
    peakwidthMAX   		= "numeric",
    ppm   				= "numeric",
    snthreshCent   		= "numeric",
    prefilterPeakNr 	= "numeric",
    prefilterIntensity	= "numeric",
    mzCenterFun 		= "character",
    integrate   		= "character",
    fitgauss   			= "logical",
    verboseColumns   	= "logical",
    noise   			= "numeric",

    ##matchedFilter arguments
    fwhm 				= "numeric",
    step   				= "numeric",
    steps   			= "numeric",
    max   				= "numeric",
    mzdiffMatched   	= "numeric",
    snthreshMatched		= "numeric",

    ##groupDensity
    minfrac   			= "numeric",
    minsamp   			= "numeric",
    bw   				= "numeric",
    mzwid   			= "numeric",	
    maxN   				= "numeric",	

    ## retcor Obiwarp
    plottype  			= "character",
    response  		 	= "numeric",
    profStep  			= "numeric",
    center  			= "numeric",
    distFunc   			= "character",
    gapInit   			= "numeric",
    gapExtend   		= "numeric",
    factorDiag   		= "numeric",
    factorGap   		= "numeric",
    localAlignment 		= "numeric",
    initPenalty   		= "numeric",

    ##retcor Peakgroups

    missing 			= "numeric",
    extra 				= "numeric",
    smooth 				= "character",
    span 				= "numeric",
    family 				= "character",
    plottypePeakgroups 	= "character",



    ##Report
    cl1   				= "character",
    cl2   				= "character",
    filename   			= "character",
    eicMax   			= "numeric",
    txtReport   		= "character",
    pdfRel   			= "character",
    pdfAbs   			= "character",
	pdfEIC       		= "character",
	pdfMZ        		= "character",
    metlinUnc  			= "numeric",
    sortpval   			= "logical",
    value  				= "character",
    eicwidth   			= "numeric",


    ##DeltaMS
    Rtwindow  			= "numeric",
    ppmw   				= "numeric",
    noiseCutoff  		= "numeric",
    monoTol   			= "logical",
    enriTol    			= "numeric",
    intChoice  			= "character",
    varEQ  				= "logical",
    alpha  				= "numeric",
    singleSample  		= "logical",
    compareOnlyDistros  = "logical",
    massOfLabeledAtom  	= "numeric",
	condition1 			= "character",
	condition2 			= "character",
	numAtom 			= "numeric",
	dpeak 				= "numeric",
	maxLT 				= "numeric",
	iRatio   			= "numeric",
    errRatio   			= "numeric",
	maxSD
	="numeric",

	##Experiment 3 specific

	uC 					= "character",
	lC 					= "character",
	uT 					= "character",
	lT 					= "character",
	class1sampNames 	= "character",
	class2sampNames 	= "character",
	sNctrl 				= "character",
	sNpert 				= "character",
	classesX13ctrl 		= "character",
	classesX13pert 		= "character",

	# Experiment 2 specific

	unLabeledS 			= "character",
	LabeledS   			= "character"

  ),


  prototype = list(

	#Additional variables

	MSFiles     		= "",
	MSfilePattern 		= "mzxml|mzxml|cdf|mzdata",
	NrFo				= NULL,
	NrFi        		= 0,
	restore 			= FALSE,
	folderNames 		= c(""),
	saveFileName 		= "",
	CorrectFilePath 	= "",
	CorrectFileName 	= "",
	WD           		= "",
	foldersOld   		= "",

    #General settings
  	higherAbundIso		= 1,
  	lowerAbundIso		= 1,
    profmethod  		= "bin",
    profparam  			= NULL,
    mslevel    			= 1,
    nSlaves  			= 1,
    TypeOfExp       	= 1,
    methodPD 			= "centWave",
    isotopeMassDiff   	= 1.003355,
    retMeth  			= "obiwarp",
	hp           		= 480, 		#Blot height EIC
	wp           		= 640,		#Blot width EIC
	iOdm        	   	= "Select two isotopes",
	manOrCalc    	 	= "Use naturally occuring ratio",
  	charge				= 1,
  	isotope1      		= "1 H",
  	isotope2      		= "2 D",


    ##centWave arguments

    mzdiffCent    		= -0.001,
    lockMassFreq   		= FALSE,
    polarity   			= NULL,
    peakwidthMIN		= 5,
    peakwidthMAX		= 20,
    ppm    				= 2.5,
    snthreshCent   		= 10,
    prefilterPeakNr 	= 3,
    prefilterIntensity 	= 5000,
    mzCenterFun  		= "wMean",
    integrate    		= "Use filtered data",
    fitgauss    		= FALSE,
    verboseColumns    	= FALSE,
    snthreshCent    	= 10,
    noise    			= 0,

    ##matchedFilter arguments

    fwhm  				= 30,
    step    			= 0.1,
    steps    			= 2,
    max    				= 5,
    mzdiffMatched    	= 0.6,
    snthreshMatched 	= 10,

    ##groupDensity
    minfrac    			= 0.5,
    minsamp    			= 1,
    bw    				= 2,
    mzwid   			= 0.015,
    maxN    			= 50,

    # retcor Obiwarp
    plottype  			= "none",
    response    		= 1,
    profStep    		= 1,
    center   			= 0,
    distFunc    		= "cor_opt",
    gapInit    			= 0,
    gapExtend    		= 0,
    factorDiag    		= 2,
    factorGap    		= 1,
    localAlignment  	= 0,
    initPenalty    		= 0,

    ##retcor Peakgroups

    missing 			= 1,
    extra 				= 1,
    smooth 				= "loess",
    span 				= 0.2,
    family 				= "gaussian",
    plottypePeakgroups 	= "none",


    ##Report
    cl1    			  	= "eg. C12",
    cl2    			  	= "eg. C13",
    filename    		= "Peaktable",
    eicMax    			= 5,
    txtReport    		= "Identified_Features",
    pdfRel    			= "Features_Relative_Intensity",
    pdfAbs    			= "Features_Absolute_Intensity",
  	pdfEIC          	= "EICs",
  	pdfMZ           	= "Mass Spectra",
    metlinUnc  			= 0.15,
    sortpval    		= TRUE,
    value   		  	= "into",
    eicwidth    		= 200,


    ##DeltaMS
    Rtwindow   			= 10,
    ppmw    			= 5,
    noiseCutoff   		= 10,
    monoTol    			= FALSE,
    enriTol     		= 0.1,
    intChoice   		= "into",
    varEQ   			= FALSE,
    alpha   			= 0.05,
    singleSample   		= FALSE,
    compareOnlyDistros 	= FALSE,
    massOfLabeledAtom  	= 12,
	condition1 			= "Control",
	condition2 			= "Perturbation",
	numAtom 			= 9,
	dpeak 				= 12,
	maxLT 				= 4,
	iRatio    			= c(1,1),
	errRatio    		= 10,
	maxSD  				= 10,

	# Experiment 3 specific

	uC 					= NULL,
	lC 					= NULL,
	uT 					= NULL,
	lT 					= NULL,
	class1sampNames 	= "",
	class2sampNames 	= "",
	sNctrl 				= "",
	sNpert 				= "",
  	classesX13ctrl 		= c("labeled", "unlabeled"),
	classesX13pert 		= c("labeled", "unlabeled"),


	# Experiment 2 specific

	unLabeledS 			= "",
	LabeledS   			= ""
  )
)


## Source functions from DeltaMSshinyappFunctions/..

source("data/DeltaMSshinyappFunctions/BtnAct.R")
source("data/DeltaMSshinyappFunctions/CheckFolders.R")
source("data/DeltaMSshinyappFunctions/iRatioFunc.R")
source("data/DeltaMSshinyappFunctions/reassignExp.R")
source("data/DeltaMSshinyappFunctions/ttWrap.R")
