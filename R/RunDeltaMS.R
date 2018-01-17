#'Start DeltaMS
#'
#'@export
#'\code{RunDeltaMS} will open DeltaMS in system's default browser.


RunDeltaMS<- function() {
  appDir <- system.file("ShinyApp", "DeltaMS", package = "DeltaMS")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing `DeltaMS`.", call. = FALSE)
  }

  shiny::runApp(appDir, display.mode = "normal", launch.browser = TRUE)
}
