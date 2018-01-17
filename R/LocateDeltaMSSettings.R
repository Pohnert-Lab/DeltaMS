#' Locate the shiny bookmarksfolder that contains saved parameter settings
#'
#'@export
#' \code{LocateDeltaMSSettings} will return the path where saved DeltaMS settings reside.
#'
#' @param number Show number of saved settings (boolean, default is FALSE)
#' @param names Show names of saved settings (boolean, default is FALSE)
#' @details \code{LocateDeltaMSSettings} optionally returns the number and the names of the present settings. If there are files of a certain setting missing, a message will be returned.

LocateDeltaMSSettings<-function(number=FALSE, names=FALSE){
pF<-file.path(.libPaths(),"DeltaMS","ShinyApp","DeltaMS","shiny_bookmarks")
message("Saved settings reside in:\n")
message(paste0(pF[which(dir.exists(file.path(.libPaths(),"DeltaMS","ShinyApp","DeltaMS","shiny_bookmarks")))],"\n"))
lenSet<-list.files(pF, pattern = ".rds", recursive = TRUE,include.dirs = FALSE,ignore.case = TRUE,full.names = TRUE)
dirsSet<-list.dirs(pF,full.names = TRUE,recursive = FALSE)
if(length(dirsSet)!=0){
if((length(lenSet)/length(dirsSet))==2){
if(number==TRUE){
message(paste(length(dirsSet),"setting(s) is/are present."))
}
if(names==TRUE){
message("The names of these settings are:\n")
basename(dirsSet)
}}else{
for (i in length(dirsSet)){
  if(!"input.rds" %in% list.files(dirsSet[i], pattern = ".rds")){
    message("The input.rds is missing in the ",basename(dirsSet[[i]])," folder.")}
  if(!"values.rds" %in% list.files(dirsSet[i], pattern = ".rds")){
    message("The values.rds is missing in the ",basename(dirsSet[[i]])," folder.")}
}
}}else{
message("There are no saved settings up to now.")
}
}