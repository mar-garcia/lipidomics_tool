library(tidyverse)
library(MetaboCoreUtils)
source("lipid_functions.R")

cmps_db <- c()

cls <- c("FA", "pHexFA", "pPentFA", "CAR", "LPA", "LPC", "LPE", "LPG", "LPI", "LPS", 
         "DGMG", "MG", "ST", "Glc-ST")
C <- seq(from = 12, to = 24, by = 1)
db <- seq(from = 0, to = 6, by = 1)
sn <-  paste(expand.grid(C, db)[,"Var1"], expand.grid(C, db)[,"Var2"], 
             sep = ":")
for(i in cls){
  cmps_db <- c(cmps_db, paste(i, sn))
}

cls <- c("SM", "Cer", "Cer;O3", "Cer;O4", "HexCer", "HexCer;O3", "HexCer;O4", "LactCer", 
         "PA", "mPA", "dmPA", "PC", "PE", "PG", "PI", "PS", "MGDG", "DGDG", "SQDG", "DGTS", "DGGA", "DG")
C <- seq(from = 12*2, to = 28*2, by = 1)
db <- seq(from = 0, to = 6*2, by = 1)
sn <-  paste(expand.grid(C, db)[,"Var1"], expand.grid(C, db)[,"Var2"], 
             sep = ":")
for(i in cls){
  cmps_db <- c(cmps_db, paste(i, sn))
}

cls <- c("acMGDG", "TG", "TG;O", "TG;O2")
C <- seq(from = 10*3, to = 24*3, by = 1)
db <- seq(from = 0, to = 6*3, by = 1)
sn <-  paste(expand.grid(C, db)[,"Var1"], expand.grid(C, db)[,"Var2"], 
             sep = ":")
for(i in cls){
  cmps_db <- c(cmps_db, paste(i, sn))
}


cmps_db <- data.frame(
  name = cmps_db,
  class = gsub("\\ .*", "", cmps_db),
  C = gsub(":.*", "", gsub(".*\\ ", "", cmps_db)),
  db = gsub(".*:", "", gsub(".*\\ ", "", cmps_db))
)
cmps_db$C <- as.numeric(cmps_db$C)
cmps_db$db <- as.numeric(cmps_db$db)
cmps_db$formula <- fml_maker(cmps_db$class, cmps_db$C, cmps_db$db)
cmps_db$mass <- calculateMass(cmps_db$formula)

save(cmps_db, file = "lipid_workspace.RData")
