### Function to open reassign menu for Experiments
### Experiment 1 isotope signature
reassignExp<-function(Nr){
  if(Nr==1 | Nr==2 | Nr==3){
  if(Nr==1){
  showModal(modalDialog(strong("Assign labels for plots"),hr(),
                                                div(style="display:inline-block",textInput("cl1", "Caption plotlabel isotope 1",value = isoD@cl1)),
                                                div(style="display:inline-block",textInput("cl2", "Caption plotlabel isotope 2",value = isoD@cl2)),

                                                footer = modalButton("Ok"),
                                                easyClose = FALSE))
    }

  ###Experiment Type 2 dataset comparison
  if(Nr==2){

  showModal(modalDialog(strong("Which sampleclass contains labeled samples?"),hr(),
                                               radioButtons("labeledFolder",label = NULL, choices = isoD@folderNames),hr(),
                                               textInput("cl1", "Caption plotlabel labeled sampleclass",value = isoD@cl2),
                                               textInput("cl2", "Caption plotlabel unlabeled sampleclass",value = isoD@cl1),

                                                footer = modalButton("Ok"),
                                                easyClose = FALSE))}


###Experiment Type 3 isotope-guided perturbation analysis
  if(Nr==3){
  showModal(modalDialog(strong("Assign sampleclasses individually"),hr(),
                                                tagList(column(6,
                                                  selectInput("lC", strong("labeled Control"),choices = c("",isoD@folderNames)),
                                                  disabled(selectInput("uC", strong("unlabeled Control"),choices = c(""))),
                                                  disabled(selectInput("lT", strong("labeled Treatment"),choices = c(""))),
                                                  disabled(selectInput("uT", strong("unlabeled Treatment"),choices = c("")))),

                                                column(6,
                                                  textInput("condition1", strong("Caption Condition 1"),value = isoD@condition1),
                                                  textInput("condition2", strong("Caption Condition 2"),value = isoD@condition2))),

                                                footer = modalButton("Ok"),
                                                easyClose = FALSE))}
  }else{
    return("Argument Nr has to be 1,2 or 3!")
  }
}



