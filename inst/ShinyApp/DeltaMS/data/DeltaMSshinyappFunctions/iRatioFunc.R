####################################################
## Ratio function to calculate isotopic ratios ###
####################################################

iRatioFunc<- function(isotope1, isotope2){

  br()
  # Check if both of the chosen isotopes are the same element (isolist.csv/isolist.RData used for this step)
  if(isoL$symbol[which(isoL$mersymbol==isotope1)]==isoL$symbol[which(isoL$mersymbol==isotope2)]) {

	if (isoL$mersymbol[which(isoL$mersymbol==isotope1)] != isoL$mersymbol[which(isoL$mersymbol==isotope2)]){

	# Determine which of the selected isotopes has the higher abundance.
    isoD@higherAbundIso<<-max(isoL$abund[which(isoL$mersymbol==isotope1)],isoL$abund[which(isoL$mersymbol==isotope2)])
	# Determine which of the selected isotopes has the lower abundance.
    isoD@lowerAbundIso<<-min(isoL$abund[which(isoL$mersymbol==isotope1)],isoL$abund[which(isoL$mersymbol==isotope2)])

	# If chosen isotope1 has the higher abundance calculate iRatio as follows and tell user the ratio
    if(isoD@higherAbundIso==isoL$abund[which(isoL$mersymbol==isotope1)]) {
      isoD@iRatio<<-c(isoD@higherAbundIso, isoD@lowerAbundIso)
      renderUI({
        tagList(
          "The ratio is", round(isoD@higherAbundIso/isoD@lowerAbundIso, digits=2), "(", strong(isotope1), ") to 1 (" , strong(isotope2),")", br(),
          isoD@higherAbundIso,"%", strong(isotope1), br(),
          isoD@lowerAbundIso,"%", strong(isotope2)
        )})

    }
    else { #If chosen isotope2 has the higher abundance calculate iRatio as follows and tell user the ratio
      isoD@iRatio<<-c(isoD@lowerAbundIso, isoD@higherAbundIso)
      renderUI({tagList(
        "The ratio is", round(isoD@higherAbundIso/isoD@lowerAbundIso, digits=2), "(",strong(isotope2),") to 1 (", strong(isotope1),")", br(),
        isoD@higherAbundIso, "%", strong(isotope2), br(),
        isoD@lowerAbundIso, "%", strong(isotope1)
      )})
    }

    }else{renderUI({
      strong("Chosen isotopes are identical")})
	  }}
  else{ # If elements do not match, no ratio will be calculated
    renderUI({
      strong("For ratio determination the same element has to be selected")
    })

  }
}

