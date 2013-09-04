library(infotheo, lib.loc="~/R/x86_64-unknown-linux-gnu-library/")
library(minet, lib.loc="~/R/x86_64-unknown-linux-gnu-library/")

# green tea
load("/protected/projects/pulmarray/Airway_Metanalysis/kyoto/aracne/data/grt_MIM.rda")  # grt.mim
grt.clr.network <- clr(grt.mim)
save(grt.clr.network, file="/protected/projects/pulmarray/Airway_Metanalysis/kyoto/aracne/data/grt_CLR.rda")
