library(infotheo, lib.loc="~/R/x86_64-unknown-linux-gnu-library/")
library(minet, lib.loc="~/R/x86_64-unknown-linux-gnu-library/")

# airway
load("/protected/projects/pulmarray/Airway_Metanalysis/kyoto/aracne/data/airway_MIM.rda")  # air.mim
air.clr.network <- clr(air.mim)
save(air.clr.network, file="/protected/projects/pulmarray/Airway_Metanalysis/kyoto/aracne/data/airway_CLR.rda")
