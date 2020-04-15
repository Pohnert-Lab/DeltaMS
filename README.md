## DeltaMS installation and base functions

To install DeltaMS copy-paste the following code:

```
if(!"devtools" %in% installed.packages()[,"Package"]){
    install.packages("devtools")}
	
library(devtools)

install_github("Pohnert-Lab/DeltaMS")
```

Afterwards load the package and check for missing packages which will be installed automatically.
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

