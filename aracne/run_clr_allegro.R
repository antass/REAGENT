library(infotheo, lib.loc="~/R/x86_64-unknown-linux-gnu-library/")
library(minet, lib.loc="~/R/x86_64-unknown-linux-gnu-library/")

# airway
load("/protected/projects/pulmarray/Airway_Metanalysis/kyoto/aracne/data/allegro16_MIM.rda")
all.clr.network <- clr(all.mim)
save(all.clr.network, file="/protected/projects/pulmarray/Airway_Metanalysis/kyoto/aracne/data/allegro_CLR.rda")
