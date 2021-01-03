###DeltaMS
library(xcms)
library(shiny)
library(shinyBS)
library(shinyjs)
library(shinythemes)
library(DescTools)
library(rChoiceDialogs)
library(multtest)
library(rmarkdown)
library(X13CMS)
source("data/SourceFile.R")
load("data/isolist.RData", envir = .GlobalEnv)
load("data/InputDefault.RData",envir = .GlobalEnv)
load("data/tooltips.RData", envir = .GlobalEnv)
isoDAlt<<-isoD

shinyApp(enableBookmarking = "server", ui = function(req){
    fluidPage(titlePanel(img(src="Penrose.jpg", width="15%"), windowTitle = "DeltaMS"),
                 shinyjs::extendShinyjs(text = "shinyjs.winprint = function() {window.print();}", functions = c("winprint")),
                 inlineCSS(".green {background: #97a7f2; font-weight: bold;}"),
                 theme = shinytheme("paper"),
                 ttWrap(),

  sidebarLayout(sidebarPanel(width = 3,
                  shinyjs::useShinyjs(),
                  shinyjs::hidden(actionButton("chDir", "Select folder", icon = icon("folder-open", lib = "glyphicon"))),
                  actionButton("aboutDeltaMS", " About DeltaMS", icon = icon("info-sign", lib = "glyphicon")),
                  p(),
                  shinyjs::hidden(actionButton("bookmBtn","Save current state", icon = icon("floppy-disk", lib = "glyphicon"))),
                  actionButton("loadB", " Load former state", icon = icon("folder-open", lib = "glyphicon")),
                  textOutput("SaveLoadText"),
                  hr(),textOutput("WD",inline=F),br(),

                  div(style="display:inline-block;vertical-align:bottom;",uiOutput("FF",inline=F)),
                  div(style="display:inline-block;vertical-align:bottom;",uiOutput("reassignBT")),hr(),


                    shinyjs::hidden(
                     radioButtons("iOdm",label = NULL, c("Select two isotopes", "Set \u0394 m manually"), inline=FALSE),
                     uiOutput("UIiOdm")),
                     shinyjs::hidden(actionButton("startAnalysis", "Start Analysis", icon = icon("play-circle", lib = "glyphicon")))
                   ),




             mainPanel(
                   tabsetPanel(id = "tabs",
                     tabPanel(strong("Select Type of Experiment"), value = "TE",fluidRow(
                       column(4, p(),
                              img(src="exp1.jpg",  width="100%"),br(),p(),
                              actionButton("exp1", "Isotope signature", icon = icon("erlenmeyer-flask", lib = "glyphicon")),p(),
                              img(src="exp3.jpg", width="207%"),br(),p(),
                              actionButton("exp3", "Isotope-guided perturbation analysis")),
                       column(4,p(),
                              img(src="exp2.jpg",width="100%"),br(),p(),
                              actionButton("exp2", "Dataset comparison"))

                     )),

                     tabPanel(strong("Report Settings"), value = "RS",fluidRow(
                       column(6,br(),p(),
                              textInput("filename", "Feature table as CSV", value = isoD@filename),
                              # textInput("txtReport", "Name of .txt report",value = isoD@txtReport),
                              textInput("pdfRel", "Bar plots (rel. Intensities as PDF)",value = isoD@pdfRel),
                              textInput("pdfAbs", "Bar plots (abs. Intensities as PDF)",value = isoD@pdfAbs)),


                       column(6,br(),p(),
                              textInput("pdfMZ", "Mass spectra as PDF",value = isoD@pdfMZ),
                              textInput("pdfEIC", "EICs as PDF",value = isoD@pdfEIC)
                              # numericInput("metlinUnc", "Metlin uncertainty", value = isoD@metlinUnc),
                              # checkboxInput("sortpval", strong("Sort results by p-value"),value = isoD@sortpval),
                              # selectInput("value", "Intensity values to be used for diffreport", c("into", "maxo", "intb"),selected = isoD@value)
                              )
                     ),hr(),
                     actionButton("ResetTabRS", "Reset inputs to default or loaded state")),

                     tabPanel(strong("Peak Detection"), value = "PD",p(),

                              radioButtons("methodPD", "Choose centWave or matchedFilter:", c("centWave", "matchedFilter")),br(),
                              uiOutput("centORmatch"),
                              hr(),
                              actionButton("ResetTabPD", "Reset inputs to default or loaded state")
                     ),

                     tabPanel(strong("Grouping"),value = "GR",p(),
                              numericInput("bw", "bw",isoD@bw),
                              numericInput("mzwid", "mzwid",isoD@mzwid),
                              numericInput("minfrac", "minfrac", isoD@minfrac),
                              numericInput("minsamp", "minsamp",isoD@minsamp),
                              numericInput("maxN", "max",isoD@maxN),hr(),
                              actionButton("ResetTabGR", "Reset inputs to default or loaded state")
                    ),
                    tabPanel(strong("Retention Time Correction"), value = "RT",p(),
                             textOutput("Exp1File1"),br(),
                             radioButtons("retMeth", "Choose method for retention time correction:", c("obiwarp", "peakgroups")),br(),
                             uiOutput("retcorOut"),hr(),
                             actionButton("ResetTabRT", "Reset inputs to default or loaded state")
                   ),
                    tabPanel(strong("DeltaMS Settings"), value = "DS",fluidRow(
                            column(6,br(),p(),
                            numericInput("Rtwindow", "Rtwindow", isoD@Rtwindow),
                            numericInput("ppmw", "ppmw", isoD@ppmw),
                            numericInput("noiseCutoff", "noiseCutoff",isoD@noiseCutoff),
                            checkboxInput("monoTol", strong("monoTol"), value = isoD@monoTol),
                            numericInput("enriTol", "enriTol", isoD@enriTol),
                            numericInput("maxLT", "maxLT", isoD@maxLT)),
                            column(6,br(),p(),
                            numericInput("dpeak", "dpeak", isoD@dpeak),
                            selectInput("intChoice", "intChoice", c("intb", "into", "maxo"), selected = isoD@intChoice),
                            checkboxInput("varEQ", strong("varEQ"), value = isoD@varEQ),
                            numericInput("alpha", "alpha", isoD@alpha),
                            numericInput("numAtom", "numAtom", isoD@numAtom),
                            checkboxInput("compareOnlyDistros", strong("compareOnlyDistros"), value = isoD@compareOnlyDistros))),
                            hr(),
                            actionButton("ResetTabDS", "Reset inputs to default or loaded state"))
                  )) #end of mainpanel

               )# end of sidebar layout
  )},



  server = function(input, output, session){
    useShinyjs()

    #About DeltaMS
    observeEvent(input$aboutDeltaMS, {
      showModal(
        modalDialog(title = "About DeltaMS",
                    HTML("Version 1.0.0 <hr>To cite the DeltaMS package in publications use:<br> 
					Baumeister, T.,Ueberschaar, N.,Schmidt-Heck, W.,Mohr, J.F.,Deicke, M.,Wichard, T.,Guthke,
					R.,Pohnert, G.<br> DeltaMS: a tool to track isotopologues in GC- and LC-MS data Metabolomics, 14:41 (2018)"),
                    footer = modalButton("Cancel"),
                    easyClose = TRUE
        ))
    })
    #Save bookmark
    observeEvent(input$bookmBtn, {
      isoD@foldersOld<<-list.dirs("shiny_bookmarks",full.names = FALSE, recursive = FALSE)
      showModal(
        modalDialog(title = "Save settings",
          textInput("saveFileName", "Enter save file name (no special characters)"),
          textOutput("FileNameAlert"),
        footer = tagList(
         disabled(actionButton("genSave", "Save")),
          modalButton("Cancel")),
        easyClose = TRUE
      ))
    })

    specialCharacter<-reactive({
      validate(
        need(!grepl("[[:punct:]]", input$saveFileName), "No special characters allowed")
      )
    })
    output$FileNameAlert <- renderText({ specialCharacter()})

    observeEvent(input$saveFileName,{
      if(!grepl("[[:punct:]]", input$saveFileName) & nzchar(input$saveFileName)){
      shinyjs::enable("genSave")
      } else {shinyjs::disable("genSave")}
      })

      observeEvent(input$genSave,{
        output$FileNameAlert<-renderText("")
        isoD@foldersOld<<-list.dirs("shiny_bookmarks",full.names = FALSE, recursive = FALSE)
        if(!any(grepl(paste0("^",input$saveFileName,"$"),isoD@foldersOld))){
          session$doBookmark()

          RecentFile<<-character(0)
          while(identical(RecentFile, character(0))){
            foldersNew<-list.dirs("shiny_bookmarks",full.names = FALSE, recursive = FALSE)
            RecentFile<<-setdiff(foldersNew, isoD@foldersOld)
          }
          f1<-file.path("shiny_bookmarks", RecentFile)
          f2<-file.path("shiny_bookmarks", isoD@saveFileName)
          AfterRenameFile<<-character(0)
          while(identical(AfterRenameFile,character(0))){
          file.rename(f1,f2)
          FoldersAfterRename<-list.dirs("shiny_bookmarks",full.names = FALSE, recursive = FALSE)
          AfterRenameFile<<-setdiff(FoldersAfterRename, foldersNew)
          }
          isoD@saveFileName<<-""
          isoD@foldersOld<<-""
          removeModal()
          output$FileNameAlert<-renderText("")

          } else {
            output$FileNameAlert<-renderText("Filename already taken")
          }
        })


    ## Bookmark exclude

    setBookmarkExclude(c("chDir",	"saveFileName","genSave","loadB","LinkGen","loadChoice","bookmBtn","loadChoice","bmLink","startAnalysis", "hashC","generateLink", "reassignSC", "ResetTabRS","ResetTabPD","ResetTabGR", "ResetTabRT", "ResetTabDS"))


    ## onBookmark
    onBookmark(function(state) {
      savedTime <- as.character(Sys.time())
      state$values$time <- savedTime
    })

    ## onBookmark Restore

    onRestore(function(state) {
    output$SaveLoadText <- renderText(paste("Restoring from state saved at:", state$values$time, "\n"))
    isoD@restore<<-TRUE

    file.rename(file.path("shiny_bookmarks", "TempFileName"),isoD@CorrectFilePath)


    # Get Input.rds file with loaded isoD settings. Otherwise isoD default various will still be used
    loadedInputSettings<-readRDS(file.path(isoD@CorrectFilePath, "input.rds"))
    LISnames<-as.list(names(loadedInputSettings))
    slotN<-slotNames(isoD)
    TempFileName<-file.path("shiny_bookmarks", "TempFileName")
    sapply(LISnames, FUN = function(i){
      if(i %in% slotN){
        slot(isoD, i)<<-loadedInputSettings[[i]]
        }
    })

    isoDAlt<<-isoD
    isoD@foldersOld<<-""

    shinyjs::alert("You restored saved settings. Select folder with MS files")
    })



    #Load bookmarked state
    observeEvent(input$loadB, {
      oldSettings<-list.dirs("shiny_bookmarks",full.names = FALSE,recursive = FALSE)
      showModal(modalDialog(

        selectInput("loadChoice", "Select saved settings file",choices = c("",oldSettings)),
        uiOutput("bmLink"),
        textOutput("MsgLoadB"),
        if(identical(oldSettings, character(0))){
          renderText("No saved Settings available")
          },
        title = "Load saved settings file",
        footer =
          actionButton("LinkGen", "Generate link"),
          modalButton("Cancel"),
        easyClose = TRUE
       ))
      })


    observeEvent(input$LinkGen,{
      isoD@CorrectFileName<<-isolate(input$loadChoice)
      isoD@CorrectFilePath<<-file.path("shiny_bookmarks",isolate(input$loadChoice))
      TempFileName<-file.path("shiny_bookmarks", "TempFileName")
      file.rename(isoD@CorrectFilePath,TempFileName)
         output$bmLink<- renderPrint(strong(a(target  = "_blank",href=link<-"/?_state_id_=TempFileName",
                                           "Load old bookmark in new Tab"),br(),
                                         strong("Close this tab after pressing link")))
      })

    observeEvent(input$loadChoice,{
      if(input$loadChoice!=""){
        if(isoD@CorrectFileName!=""){
        if(input$loadChoice!=isoD@CorrectFileName){
          file.rename(file.path("shiny_bookmarks", "TempFileName"),file.path("shiny_bookmarks",isoD@CorrectFileName))
          output$bmLink<-renderText("")
      }}}
    })

    #Reset DeltaMS
    # observeEvent(input$refreshBtn,{
    #  showModal(modalDialog(
    #    "Do you really want to reset the DeltaMS settings?",
    #     footer = tagList(
    #       modalButton("Cancel"),
    #       actionButton("resetYes","Yes, reset DeltaMS")),
    #    easyClose = TRUE
    #   ))
    # observeEvent(input$resetYes, {
    #   js$refresh()
    # })
    # })

  #   ###ObserveEvent to save inputs into isoD object
  # isoD<<-InputDefaults() # In case of destroying default DeltaMS default .RData, remove # and create new isoD
  slnames<-slotNames(isoD)
  lapply(slnames, function(i) {
  observeEvent(input[[i]], {
    if(is.numeric(slot(isoD, i))){
      if(is.numeric(input[[i]])){
        slot(isoD, i)<<-input[[i]]
      }}else{
      slot(isoD, i)<<-input[[i]]}
  })})

  ## Tab update Buttons
  observeEvent(input$ResetTabRS,{
    updateTextInput(session, "filename",  value = isoDAlt@filename)
    updateTextInput(session,"txtReport", value = isoDAlt@txtReport)
    updateTextInput(session,"pdfRel", value = isoDAlt@pdfRel)
    updateTextInput(session,"pdfAbs", value = isoDAlt@pdfAbs)
    updateTextInput(session,"pdfMZ", value = isoDAlt@pdfMZ)
    updateTextInput(session,"pdfEIC", "PDF EICs",value = isoDAlt@pdfEIC)
    updateNumericInput(session, "metlinUnc",  value = isoDAlt@metlinUnc)
    updateCheckboxInput(session, "sortpval", value = isoDAlt@sortpval)
    updateSelectInput(session, "value", selected = isoDAlt@value)
  })

  observeEvent(input$ResetTabPD,{
    if(isoD@methodPD=="centWave"){
      updateNumericInput(session,"ppm", value = isoDAlt@ppm)
      updateNumericInput(session,"peakwidthMIN",value = isoDAlt@peakwidthMIN)
      updateNumericInput(session,"peakwidthMAX",value = isoDAlt@peakwidthMAX)
      updateNumericInput(session,"snthreshCent",value = isoDAlt@snthreshCent)
      updateNumericInput(session,"prefilterPeakNr",value = isoDAlt@prefilterPeakNr)
      updateNumericInput(session,"prefilterIntensity",value = isoDAlt@prefilterIntensity)
      updateNumericInput(session,"mzdiffCent",value = isoDAlt@mzdiffCent)
      updateNumericInput(session,"noise",value = isoDAlt@noise)
      updateSelectInput(session, "polarity",selected = isoDAlt@polarity)
      updateSelectInput(session, "mzCenterFun", selected = isoDAlt@mzCenterFun)
      updateSelectInput(session, "integrate", selected = isoDAlt@integrate)
      updateCheckboxInput(session, "fitgauss", value = isoDAlt@fitgauss)
      updateCheckboxInput(session, "verboseColumns", value = isoDAlt@verboseColumns)
    }else{
      updateNumericInput(session, "fwhm", value = isoDAlt@fwhm)
      updateNumericInput(session, "mzdiffMatched", value = isoDAlt@mzdiffMatched)
      updateNumericInput(session, "step", value = isoDAlt@step)
      updateNumericInput(session, "steps", value = isoDAlt@steps)
      updateNumericInput(session, "max", value = isoDAlt@max)
      updateNumericInput(session, "snthreshMatched", value = isoDAlt@snthreshMatched)
    }
  })

  observeEvent(input$ResetTabGR,{
    updateNumericInput(session, "bw", value = isoDAlt@bw)
    updateNumericInput(session, "mzwid", value = isoDAlt@mzwid)
    updateNumericInput(session, "minfrac", value = isoDAlt@minfrac)
    updateNumericInput(session, "minsamp", value = isoDAlt@minsamp)
    updateNumericInput(session, "maxN", value = isoDAlt@maxN)
  })

  observeEvent(input$ResetTabRT,{
    if(isoD@retMeth == "obiwarp"){
      updateSelectInput(session, "plottype", selected = isoDAlt@plottype)
      updateNumericInput(session, "response", value = isoDAlt@response)
      updateNumericInput(session, "profStep", value = isoDAlt@profStep)
      updateNumericInput(session, "center", value = isoDAlt@center)
      updateSelectInput(session, "distFunc", selected = isoDAlt@distFunc)
      updateNumericInput(session, "gapExtend", value = isoDAlt@gapExtend)
      updateNumericInput(session, "gapInit", value = isoDAlt@gapInit)
      updateNumericInput(session, "factorDiag", value = isoDAlt@factorDiag)
      updateNumericInput(session, "factorGap", value = isoDAlt@factorGap)
      updateNumericInput(session, "initPenalty", value = isoDAlt@initPenalty)
      updateNumericInput(session, "localAlignment", value = isoDAlt@localAlignment)
    }else{
      updateSelectInput(session, "smooth", selected = isoDAlt@smooth)
      updateNumericInput(session, "span", value = isoDAlt@span)
      updateSelectInput(session, "family", selected = isoDAlt@family)
      updateNumericInput(session, "extra", value = isoDAlt@extra)
      updateSelectInput(session, "plottypePeakgroups", selected = isoDAlt@plottypePeakgroups)
    }
  })

  observeEvent(input$ResetTabDS,{
    updateNumericInput(session, "Rtwindow", value = isoDAlt@Rtwindow)
    updateNumericInput(session, "ppmw", value = isoDAlt@ppmw)
    updateNumericInput(session, "noiseCutoff", value = isoDAlt@noiseCutoff)
    updateCheckboxInput(session, "monoTol", value = isoDAlt@monoTol)
    updateNumericInput(session, "enriTol", value = isoDAlt@enriTol)
    updateSelectInput(session, "intChoice", selected = isoDAlt@intChoice)
    updateCheckboxInput(session, "varEQ", value = isoDAlt@varEQ)
    updateNumericInput(session, "alpha", value = isoDAlt@alpha)
    updateCheckboxInput(session, "compareOnlyDistros", value = isoDAlt@compareOnlyDistros)
  })

  ###


    ##Experiment select buttons
    expL<-c("exp1", "exp2", "exp3")
    lapply(expL, function(i) {
      observeEvent(input[[i]], {
        list<-list("saveB","bookmBtn","iOdm","UIiOdm","chDir", "WD", "FF", "nCoreOut","refreshBtn","loadB")
        lapply(list, shinyjs::show,anim=FALSE)})
      })

    observeEvent(input$exp1, {
      isoD@TypeofExp<<-1
      shinyjs::enable("isot1")
      shinyjs::enable("isot2")
      shinyjs::enable("manOrCalc")
      shinyjs::enable("errRatio")
      shinyjs::enable("maxSD")
      shinyjs::addClass("exp1","green")
      shinyjs::removeClass("exp2", "green")
      shinyjs::removeClass("exp3", "green")
      if(!is.null(isoD@NrFo)){
        if(isoD@NrFo!=0 & isoD@restore == FALSE| isoD@NrFi == 0){
          BtnAct("startAnalysis", "red", "exp")
          BtnAct("reassignSC", "red")
          }
        if(isoD@NrFo==0 & isoD@NrFi == 1 & isoD@restore == FALSE){
          BtnAct("startAnalysis")
          BtnAct("reassignSC", "no")
        }
        }
      })

    observeEvent(input$exp2, {
      isoD@TypeofExp<<-2
      shinyjs::enable("isot1")
      shinyjs::enable("isot2")
      shinyjs::enable("manOrCalc")
      shinyjs::enable("errRatio")
      shinyjs::enable("maxSD")
      shinyjs::addClass("exp2","green")
      shinyjs::removeClass("exp1", "green")
      shinyjs::removeClass("exp3", "green")
      if(!is.null(isoD@NrFo)){
        if(isoD@NrFo!=2 & isoD@restore == FALSE| isoD@NrFi == 0){
          BtnAct("reassignSC", "red")
          BtnAct("startAnalysis", "red", "exp")
          }
        if(isoD@NrFo==2 & isoD@NrFi != 0 & isoD@restore == FALSE){
          BtnAct("startAnalysis")
          BtnAct("reassignSC", "no")
          }
        }
      })

    observeEvent(input$exp3, {
      isoD@TypeofExp<<-3
      shinyjs::disable("manOrCalc")
      shinyjs::disable("isot1")
      shinyjs::disable("isot2")
      shinyjs::disable("errRatio")
      shinyjs::disable("maxSD")
      shinyjs::addClass("exp3","green")
      shinyjs::removeClass("exp1", "green")
      shinyjs::removeClass("exp2", "green")
      if("X13CMS" %in% installed.packages()[,"Package"]){
      if(!is.null(isoD@NrFo)){
        if(isoD@NrFo!=4 & isoD@restore == FALSE| isoD@NrFi== 0){
          BtnAct("reassignSC", "red")
          BtnAct("startAnalysis", "red", "exp")
          }
        if(isoD@NrFo==4 & isoD@NrFi != 0 & isoD@restore == FALSE){
          BtnAct("startAnalysis")
          BtnAct("reassignSC", "no")
          }
      }}else{
        showModal(
          modalDialog(title = "X13CMS Package is missing",
                      tags$a(target  = "_blank",href="http://pattilab.wustl.edu/software/x13cms/x13cms.php", "Download X13CMS package here"),
                      footer = modalButton("Cancel"),
                      easyClose = TRUE
          ))
        BtnAct("startAnalysis", "red")
        BtnAct("reassignSC", "no")
        }
      })


    ## Choose directory (experiment dependent)
	
	WD.react<-reactiveVal(value = getwd()) # Reactive value to change initial folder to last chosen folder.
	
    observeEvent(input$chDir, {
      if(canUseJava()){
      isoD@WD<<-jchoose.dir(caption="Select path to MS-files containing folders", default = WD.react())
      }else{
	  isoD@WD<<-tcltk::tk_choose.dir(caption="Select path to MS-files containing folders", default = WD.react())
	  }
      if (!identical(isoD@WD, character(0))) {
        WD.react(isoD@WD)
        normWDD<-normalizePath(isoD@WD)
        if(isoD@TypeofExp!=1){
        isoD@MSFiles<<-(list.files(file.path(isoD@WD,names(CheckFolders(isoD@WD,isoD@MSfilePattern))), full.names = TRUE,pattern = isoD@MSfilePattern,ignore.case = TRUE))
        }else{isoD@MSFiles<<-list.files(isoD@WD,recursive = FALSE, full.names = TRUE, pattern = isoD@MSfilePattern,ignore.case = TRUE)}
        isoD@NrFi<<-length(isoD@MSFiles)
        if(isoD@NrFi==1){
          output$Exp1File1<-renderText("Retention time correction will not be used due to lack of replicates (one MS file selected)")
        }
        isoD@NrFo<<-length(CheckFolders(isoD@WD, isoD@MSfilePattern))
        isoD@folderNames<<-names(CheckFolders(isoD@WD, isoD@MSfilePattern))
          if (isoD@TypeofExp==1) {
          if(isoD@NrFo==0){
              if (isoD@NrFi==0) {BtnAct("startAnalysis", "red", "format")
              }else{
              BtnAct("startAnalysis")

              output$WD<-renderText(normWDD)

              output$FF<-renderText(paste(strong(isoD@NrFi),"MS-files detected"))

              reassignExp(Nr = 1)

              output$reassignBT<-renderUI({
                actionButton("reassignSC", "Reassign")
              })
              }


          }else{
                BtnAct("startAnalysis", "red","folder")
                isoD@NrFi<<-0
                isoD@NrFo<<-0
            }

        }else{
          if(isoD@TypeofExp==2){
          if(isoD@NrFo==2){
              if (isoD@NrFi==0) {BtnAct("startAnalysis", "red", "format")
              }else{
              BtnAct("startAnalysis")

              reassignExp(Nr = 2)

              output$WD<-renderText(normWDD)

              output$FF<-renderUI({
              fileNumber<-paste(strong(isoD@NrFi),"MS-files detected")
              fileNrF1<-length(list.files(file.path(isoD@WD,input$labeledFolder),pattern=isoD@MSfilePattern,ignore.case = TRUE))
              fileNrF2<-length(list.files(file.path(isoD@WD,isoD@folderNames[!isoD@folderNames %in%  input$labeledFolder]),pattern=isoD@MSfilePattern,ignore.case = TRUE))
              folderLabeled<-paste(strong(fileNrF1),"Labeled samples in sampleclass:", strong(isoD@LabeledS<<-as.character(input$labeledFolder)))
              folderUnlabeled<-paste(strong(fileNrF2), "Unlabeled samples in sampleclass:", strong(isoD@unLabeledS<<-isoD@folderNames[which(input$labeledFolder!=isoD@folderNames)]))
              HTML(paste (fileNumber, folderLabeled, folderUnlabeled, sep = '<br/>' ))
              })

              output$reassignBT<-renderUI({
                actionButton("reassignSC", "Reassign")

              })

              }}else{
                BtnAct("startAnalysis", "red","folder",2)
                isoD@NrFi<<-0
                isoD@NrFo<<-0
                }


        }else{
          if(isoD@NrFo==4){
          if(isoD@NrFi==0) {BtnAct("startAnalysis", "red", "format")
                }else{
                BtnAct("startAnalysis")

                  reassignExp(Nr = 3)


                    observeEvent(input$lC, {
                    uC<-c("",c(isoD@folderNames[-which(isoD@folderNames==input$lC)]))
                    updateSelectInput(session, "uC", choices = uC)
                    shinyjs::toggleState("uC", condition = input$lC!="")
                    })

                    observeEvent(input$uC, {
                    lT<-c("",c(isoD@folderNames[-c(which(isoD@folderNames==input$lC),which(isoD@folderNames==input$uC))]))
                    updateSelectInput(session, "lT", choices = lT)
                    shinyjs::toggleState("lT", condition = input$uC!="")
                    })

                    observeEvent(input$lT, {
                    uT<-c("",c(isoD@folderNames[-c(which(isoD@folderNames==input$lC),which(isoD@folderNames==input$uC),which(isoD@folderNames==input$lT))]))
                    updateSelectInput(session, "uT", choices = uT)
                    shinyjs::toggleState("uT", condition = input$lT!="")
                    })

                    observeEvent(input$uT, {
                      if(input$uC!="" & input$lC!="" & input$uT!="" & input$lT!="") {
                        isoD@classesX13ctrl<<-c(rep(input$uC,length(isoD@class1sampNames<<-list.files(file.path(isoD@WD,input$uC), pattern=isoD@MSfilePattern, full.names = FALSE, ignore.case = TRUE))),rep(input$lC,length(list.files(file.path(isoD@WD,input$lC),pattern=isoD@MSfilePattern, full.names = FALSE, ignore.case = TRUE))))
                        isoD@classesX13pert<<-c(rep(input$uT,length(isoD@class2sampNames<<-list.files(file.path(isoD@WD,input$uT), pattern=isoD@MSfilePattern, full.names = FALSE, ignore.case = TRUE))),rep(input$lT,length(list.files(file.path(isoD@WD,input$lT),pattern=isoD@MSfilePattern, full.names = FALSE, ignore.case = TRUE))))
                      }
                    })

                  output$WD<-renderText(normWDD)


                  output$FF<-renderUI({
                    fileNumber<-paste(isoD@NrFi,"MS-files detected")
                    fileNrF1<-length(lC<-list.files(file.path(isoD@WD,input$lC),pattern=isoD@MSfilePattern, full.names = FALSE, ignore.case = TRUE))
                    fileNrF2<-length(uC<-list.files(file.path(isoD@WD,input$uC),pattern=isoD@MSfilePattern, full.names = FALSE, ignore.case = TRUE))
                    fileNrF3<-length(lT<-list.files(file.path(isoD@WD,input$lT),pattern=isoD@MSfilePattern, full.names = FALSE, ignore.case = TRUE))
                    fileNrF4<-length(uT<-list.files(file.path(isoD@WD,input$uT),pattern=isoD@MSfilePattern, full.names = FALSE, ignore.case = TRUE))

                    isoD@sNctrl<<-c(uC,lC)
                    isoD@sNpert<<-c(uT,lT)

                    folderlC<-paste(strong(fileNrF1),"Labeled control samples in sampleclass:", strong(input$lC))
                    folderuC<-paste(strong(fileNrF2),"Unlabeled control samples in sampleclass:", strong(input$uC))
                    folderlT<-paste(strong(fileNrF3),"Labeled treatment samples in sampleclass:", strong(input$lT))
                    folderuT<-paste(strong(fileNrF4),"Unlabeled treatment samples in sampleclass:", strong(input$uT))

                    HTML(paste(fileNumber, folderlC,folderuC,folderlT,folderuT, sep = '<br/>' ))})

                  output$reassignBT<-renderUI({
                    actionButton("reassignSC", "Reassign")
                    })
                  }
          }else{
                  BtnAct("startAnalysis", "red","folder",4)
                  isoD@NrFi<<-0
                  isoD@NrFo<<-0
            }

      }}}})


    ### Reassign Button

    observeEvent(input$reassignSC,{
      if(isoD@TypeofExp==1){
        reassignExp(Nr = 1)}
      if(isoD@TypeofExp==2){
        reassignExp(Nr = 2)}
      if(isoD@TypeofExp==3){
        reassignExp(Nr = 3)}
    })

    ###Tab Peakdetection centWave
    output$centORmatch <- renderUI({

      if (input$methodPD=="centWave"){fluidRow(ttWrap(),
        column(6, tagList(
                numericInput("ppm", "ppm",isoD@ppm),
                div(style="display:inline-block;vertical-align:top; width: 150px",numericInput("peakwidthMIN", label="peakwidth (min)", isoD@peakwidthMIN)),
                div(style="display:inline-block;vertical-align:top; width: 150px",numericInput("peakwidthMAX", label="peakwidth (max)", isoD@peakwidthMAX)),
                numericInput("snthreshCent", "snthresh",isoD@snthreshCent),
                div(style="display:inline-block;vertical-align:top; width: 150px",numericInput("prefilterPeakNr", label="prefilter (k peaks)", isoD@prefilterPeakNr)),
                div(style="display:inline-block;vertical-align:top; width: 150px",numericInput("prefilterIntensity", label="prefilter (intensity I)", isoD@prefilterIntensity)),
                numericInput("mzdiffCent", "mzdiff",isoD@mzdiffCent),
                numericInput("noise", "noise",isoD@noise))),
        column(6, tagList(
                selectInput("polarity", "polarity", c("positive","negative","NULL"), selected = isoD@polarity),
                selectInput("mzCenterFun", "mzCenterFun", c("wMean", "mean", "apex", "wMeanApex3", "meanApex3"),selected = isoD@mzCenterFun),
                selectInput("integrate", "integrate", c("Use filtered data","Use real data (prone to noise)"), selected = isoD@integrate),
                checkboxInput("fitgauss", strong("fitgauss"), isoD@fitgauss),
                checkboxInput("verboseColumns", strong("verbose.columns"), isoD@verboseColumns))))

    ###Tab Peakdetection matchedFilter
      }else{ tagList(ttWrap(),
        numericInput("fwhm", "fwhm",isoD@fwhm),
        numericInput("mzdiffMatched","mzdiffMatched", isoD@mzdiffMatched),
        numericInput("step", "step",isoD@step),
        numericInput("steps", "steps", isoD@steps),
        numericInput("max", "max", isoD@max),
        numericInput("snthreshMatched", "snthresh",isoD@snthreshMatched))
      }
    })


    ###Tab Retcor obiwarp
   output$retcorOut <-renderUI({

     if (input$retMeth=="obiwarp"){fluidRow(ttWrap(),
       column(6, tagList(
                selectInput("plottype", "plottype", c("deviation", "none"), isoD@plottype),
                numericInput("response", "response", isoD@response),
                numericInput("profStep", "profStep", isoD@profStep),
                numericInput("center", "center",value = isoD@center),
                selectInput("distFunc", "distFunc", c("cor_opt", "cor", "cov", "prd", "euc"), isoD@distFunc))),
      column(6, tagList(
                numericInput("gapExtend", "gapExtend", isoD@gapExtend),
                numericInput("gapInit", "gapInit", isoD@gapInit),
                numericInput("factorDiag", "factorDiag", isoD@factorDiag),
                numericInput("factorGap", "factorGap", isoD@factorGap),
                numericInput("initPenalty", "initPenalty", isoD@initPenalty),
                numericInput("localAlignment", "localAlignment", isoD@localAlignment)
      ))
      )
     }else{  fluidRow(
       ttWrap(),
       column(6, tagList(
                selectInput("smooth", "smooth", c("loess", "linear"), isoD@smooth),
                numericInput("span", "span", isoD@span),
                selectInput("family", "family", c("gaussian", "symmetric"), isoD@family)
                )),
       column(6,tagList(
                numericInput("missing", "missing", isoD@missing),
                numericInput("extra", "extra", isoD@extra),
                selectInput("plottypePeakgroups", "plottype", c("deviation", "mdevden", "none"),isoD@plottypePeakgroups)))
       )}
     })

   observeEvent(input$smooth, {
     if (input$smooth=="loess") {
       shinyjs::show(id = "span")
     }else{
       shinyjs::hide(id ="span")
     }
   })


   ###Sidepanel choices

    output$UIiOdm <- renderUI({
if(isoD@TypeofExp!=3){
  tagList(if (input$iOdm=="Select two isotopes") {
    tagList(selectInput("isotope1", "Light isotope", isoL$mersymbol, selected = isoD@isotope1, width="70%"),
            selectInput("isotope2", "Heavy isotope", isoL$mersymbol, selected = isoD@isotope2,width="70%"),

            output$delIso1Iso2<- renderUI({
              tagList(
                "The mass difference is",
                strong(round(isoD@isotopeMassDiff<<-abs(isoL$mass[which(isoL$mersymbol==input$isotope2)]-isoL$mass[which(isoL$mersymbol==input$isotope1)]), digits = 5)),
                strong(" u."),
                p(),
                iRatioFunc(input$isotope1, input$isotope2),p(),
                numericInput("charge", "Select charge z of analyzed compounds", width="70%", value =  isoD@charge),
                radioButtons("manOrCalc", label=NULL, c("Use naturally occuring ratio","Set isotope ratio manually"), selected = isoD@manOrCalc),
                renderUI({if (input$manOrCalc=="Set isotope ratio manually")
                {
                  tagList(
                    div(style="display:inline-block",numericInput("isot1", label=input$isotope1, 1, width="65%")),
                    div(style="display:inline-block",numericInput("isot2", label=input$isotope2, 1, width="65%")),ttWrap())
                }
                }),
                div(style="display:inline-block",numericInput("errRatio","Deviation of the intensity ratio (in %)", isoD@errRatio, width="65%")),
                div(style="display:inline-block",numericInput("maxSD","SD of replicate intensities (in %)", isoD@maxSD, width="65%")),
                ttWrap())
            }))

  }else{
    tagList(numericInput("isotopeMassDiff", "\u0394 m", value=1.006277, width="70%"),br(),
            numericInput("charge", "Select charge z of analyzed compounds", width="70%", value = isoD@charge),br(),
            numericInput("massOfLabeledAtom", "Exact mass of labeled Atom",2.014102, width="70%"),br(),
            div(style="display:inline-block",strong("Isotope ratio:")),br(),
            div(style="display:inline-block",numericInput("isot1", label=HTML(paste ("Isotope 1","(unlabeled sample)", sep = '<br/>' )), 1, width="65%")),
            div(style="display:inline-block",numericInput("isot2", label=HTML(paste ("Isotope 2","(labeled sample)", sep = '<br/>' )), 1, width="65%")),br(),

            div(style="display:inline-block",numericInput("errRatio","Deviation of the intensity ratio (in %)", isoD@errRatio, width="65%")),
            div(style="display:inline-block",numericInput("maxSD","SD of replicate intensities (in %)", isoD@maxSD, width="65%")),
            ttWrap())

  })


     } else{
        tagList(if (input$iOdm=="Select two isotopes") {
          tagList(selectInput("isotope1", "Unlabeled sample", isoL$mersymbol, selected = isoD@isotope1, width="70%"),
                  selectInput("isotope2", "Labeled Sample", isoL$mersymbol, selected = isoD@isotope2,width="70%"),

                  output$delIso1Iso2<- renderUI({
                    tagList(
                      "The mass difference is",
                      strong(round(isoD@isotopeMassDiff<<-abs(isoL$mass[which(isoL$mersymbol==input$isotope2)]-isoL$mass[which(isoL$mersymbol==input$isotope1)]), digits = 5)),
                      strong(" u."),
                      p(),
                      iRatioFunc(input$isotope1, input$isotope2),p(),
                      numericInput("charge", "Select charge z of analyzed compounds", width="70%", value =  isoD@charge),
                      disabled(radioButtons("manOrCalc", label=NULL, c("Use naturally occuring ratio","Set isotope ratio manually"), selected = isoD@manOrCalc)),
                      renderUI({
                        if (input$manOrCalc=="Set isotope ratio manually"){
                          tagList(
                            div(style="display:inline-block",numericInput("isot1", label=input$isotope1, 1, width="65%")),
                            div(style="display:inline-block",numericInput("isot2", label=input$isotope2, 1, width="65%")),ttWrap())
                        }}),
                      renderUI({
                        tagList(
                          div(style="display:inline-block",disabled(numericInput("errRatio","Deviation of the intensity ratio (in %)", isoD@errRatio, width="50%"))),
                          div(style="display:inline-block",disabled(numericInput("maxSD","Max. RSD of replicate intensities (in %)", isoD@maxSD, width="50%"))),
                          ttWrap())}))
                  }))

        }else{
          tagList(numericInput("isotopeMassDiff", "\u0394 m", value=1.006277, width="70%"),br(),
                  numericInput("charge", "Select charge z of analyzed compounds", width="70%", value = isoD@charge),br(),
                  numericInput("massOfLabeledAtom", "Exact mass of labeled Atom",2.014102, width="70%"),br(),
                  div(style="display:inline-block",strong("Isotope ratio:")),br(),
                  div(style="display:inline-block",disabled(numericInput("isot1", label=HTML(paste ("Isotope 1","(unlabeled sample)", sep = '<br/>' )), 1, width="65%"))),
                  div(style="display:inline-block",disabled(numericInput("isot2", label=HTML(paste ("Isotope 2","(labeled sample)", sep = '<br/>' )), 1, width="65%"))),br(),

                  div(style="display:inline-block",disabled(numericInput("errRatio","Deviation of the intensity ratio (in %)", isoD@errRatio, width="65%"))),
                  div(style="display:inline-block",disabled(numericInput("maxSD","SD of replicate intensities (in %)", isoD@maxSD, width="65%"))),
                  ttWrap())

        })
      }
      })

    #observeEvent for labeled isotope mass
    observeEvent(input$isotope1, {
      isoD@massOfLabeledAtom<<- (isoL$mass[which(isoL$mersymbol==input$isotope1)]-1)
    })

    # if manual setting of isotope ratios
    observeEvent(input$isot1,{
      if(is.numeric(input$isot1)){
        isoD@iRatio[1]<<-input$isot1
      }
    })

    observeEvent(input$isot2,{
      if(is.numeric(input$isot2)){
        isoD@iRatio[2]<<-input$isot2
      }
    })

    observeEvent(input$manOrCalc, {
      if(input$manOrCalc=="Use naturally occuring ratio") {
        isoD@iRatio<<-c(isoD@higherAbundIso, isoD@lowerAbundIso)
        }
      if(input$manOrCalc=="Set isotope ratio manually"){
        if(is.numeric(c(input$isot1,input$isot2))){
        isoD@iRatio<<-c(input$isot1,input$isot2)
        }
      }
    })

    observeEvent(input$iOdm,{
      if(input$iOdm=="Set \u0394 m manually"){
        if(is.numeric(c(input$isot1,input$isot2))){
          isoD@iRatio<<-c(input$isot1,input$isot2)
        }
      }
    })


    ##exp2 labeled&unlabeled samples test
    observeEvent(input$labeledFolder, {
    isoD@LabeledS<<-isoD@folderNames[which(isoD@folderNames==input$labeledFolder)]
    isoD@unLabeledS<<-isoD@folderNames[which(isoD@folderNames!=input$labeledFolder)]
    })


    ## numAtom and dpeak has to be >=2, therefore input will be set back to 2 if input is < 2
    observeEvent(input$numAtom, {
    if(input$numAtom<2){
    isoD@numAtom<<-2
    updateNumericInput(session, "numAtom", value = 2)}
    })

    observeEvent(input$dpeak, {
      if(input$dpeak<2){
        isoD@dpeak<<-2
        updateNumericInput(session, "dpeak", value = 2)}
    })

###startAnalysis button
observeEvent(input$startAnalysis,{
  cat("\014")
  
  shinyjs::disable(id = "startAnalysis") # Disable Start analysis button until analysis is finished to prevent multiple activations by chance

  # iRatio depending on manually set or natural isotope ratio chosen
  if(isoD@manOrCalc != "Use naturally occuring ratio"){
    isoD@iRatio<<-round((isoD@iRatio/sum(isoD@iRatio))*100,digits = 2)
  }

  # Results folder
  resultsDir<<-file.path(isoD@WD, paste0("DeltaMS_results_",format(Sys.time(), "%m%d%y-%H%M%p")))
  if(!dir.exists(resultsDir)){
    dir.create(resultsDir)
  }else{
    dir.create(paste0(resultsDir, "_1"))
  }

  # Render Analysis_settings html file
  rmarkdown::render("data/Markdown/MainScript.Rmd",
  output_file = "Analysis_settings.html", output_dir = resultsDir,
  envir = new.env(parent = globalenv()))


if(isoD@TypeofExp==1){

  source("data/functionsDeltaMS.R")
  source("data/exp1.R")
}

  if(isoD@TypeofExp==2){
    source("data/functionsDeltaMS.R")
    source("data/exp2.R")
  }


if(isoD@TypeofExp==3){
  # source("data/functionsv6.R")
  source("data/exp3.R")
}
  
  shinyjs::enable(id = "startAnalysis")  # enable start analysis button again after analysis has finished
  
})

##End of App
})











