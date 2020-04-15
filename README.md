## DeltaMS installation and base functions

To install DeltaMS copy-paste the following code:

```
if(!"devtools" %in% installed.packages()[,"Package"]){
    install.packages("devtools")}
	
library(devtools)

install_github("Pohnert-Lab/DeltaMS")
```

Afterwards load the package and check for missing packages which will be installed automatically.<br>
If a package cannot be installed, try to manually install it via e.g. CRAN (install.packages("PackageName").

```
library(DeltaMS)
MissingPackagesDeltaMS()
```
Execute MissingPackagesDeltaMS() again to confirm that all packages seem present.

Start DeltaMS

```
RunDeltaMS()
```

The location of the saved DeltaMS settings can be found with this command

```
LocateDeltaMSSettings()
```
## Notes
Last tested versions: <br>
R version 3.6.1 <br>
RStudio 1.2.5042 (older versions may have trouble with pandoc.exe (conversion error 5) that give the html analysis file output)

RJava runs into troubles if installed Java version is not matching RStudio in terms of 32/64bit <br>
Java 64bit has to be downloaded manually on <br>
https://www.java.com/de/download/manual.jsp <br>
