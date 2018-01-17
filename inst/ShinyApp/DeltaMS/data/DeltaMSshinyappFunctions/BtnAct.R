######################################
# Folder Button colour and message ##
######################################
#This function is responsible for changing the activity and colour of the analysis button

BtnAct <- function(btn,colour ="green", msg="noMessage", nrFo = 1) {
  if (colour=="green") {tagList(
	shinyjs::show(btn),
    shinyjs::addClass(btn, "green"),
    shinyjs::enable(btn),
	  isoD@restore<<-FALSE)}
  if(colour=="red"){tagList(
    shinyjs::removeClass(btn, "green"),
    shinyjs::disable(btn))
    }
  else{
    shinyjs::enable(btn)}

  if(msg=="format"){
    if(isoD@TypeofExp==1){
      shinyjs::alert("Only one converted raw-file allowed for this analysis")
    }else{
    shinyjs::alert("No files in suitable raw data format detected. Please select another folder")
    }}

  if(msg=="folder"){
    shinyjs::alert(paste("DeltaMS expects", nrFo, "folder(s) with MS files"))
    }

	if(msg=="exp"){
      shinyjs::alert("Change folder or type of analysis")}
}
