###Runs all scripts required for the Reef Life Explorer
## EClausius November 2020 

try(
  source("R/DFs for RLE.R"), 
  TRUE)

try(source("R/Country IDW for RLE.R"), 
    TRUE)

try(source("R/Location indicator GAMs for RLE.R"), 
    TRUE)

try(source("R/Country IDW for RLE.R"), 
    TRUE)

try(source("R/Regional Reef Health for RLE.R"), 
    TRUE)
