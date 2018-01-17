############################
# tooltips wrapper ######
############################
# This function uses the tooltips.csv to render all the tooltips.
ttWrap<- function(){
  t<-tooltips
  lapply(seq(nrow(t)), function (i) {bsTooltip(id=as.character(t$id[i]),
                                               title = as.character(t$text[i]),
                                               placement = as.character(t$position[i]),
                                               trigger = as.character(t$action[i]))}
  )}
