library(infotheo, lib.loc="~/R/x86_64-unknown-linux-gnu-library/")
library(minet, lib.loc="~/R/x86_64-unknown-linux-gnu-library/")

# green tea
load("/protected/projects/pulmarray/Airway_Metanalysis/kyoto/aracne/data/lam_MIM.rda")  # lam.mim
lam.clr.network <- clr(lam.mim)
save(lam.clr.network, file="/protected/projects/pulmarray/Airway_Metanalysis/kyoto/aracne/data/lam_CLR.rda")
