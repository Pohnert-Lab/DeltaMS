#'Install DeltaMS dependent packages
#'
#'@export
#'\code{MissingPackagesDeltaMS} will install DeltaMS dependent packages except X13CMS which is required for the "Isotopologue-guided perturbation analysis"
#'@details A list of necessary packages will be checked against your installed packages. Missing packages will be installed.The X13CMS package has to be installed from an external \href{http://pattilab.wustl.edu/software/x13cms/x13cms.php}{source}.

MissingPackagesDeltaMS<-function(){
  if(!"xcms" %in% installed.packages()[,"Package"]){
    source("https://bioconductor.org/biocLite.R")
    biocLite("xcms")}
  if(!"X13CMS" %in% installed.packages()[,"Package"]){
    message("X13CMS package is missing. \nTo install X13CMS visit http://pattilab.wustl.edu/software/x13cms/x13cms.php")}
  list.of.packages <- c("shiny","shinyjs","shinyBS","shinythemes","DescTools", "rChoiceDialogs", "multtest", "rmarkdown","V8","rJava")
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)){
    install.packages(new.packages)}
  if("X13CMS" %in% installed.packages()[,"Package"]&&identical(new.packages, character(0))){
    message("All necessary packages seem to be present.")}
  }