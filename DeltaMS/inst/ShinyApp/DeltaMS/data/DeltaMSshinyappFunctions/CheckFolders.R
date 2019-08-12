############################################################
## CheckFolder function to count folders with MS files##
############################################################
CheckFolders<- function (WD, pattern) {
  dirs<-list.dirs(WD,recursive = FALSE, full.names = TRUE)

  l<-lapply(dirs, FUN = list.files, pattern = pattern, ignore.case = TRUE)
  names(l)<-list.dirs(WD,recursive = FALSE, full.names = FALSE)
  l[l!="character(0)"]
}
