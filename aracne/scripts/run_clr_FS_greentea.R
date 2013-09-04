library(infotheo, lib.loc="~/R/x86_64-unknown-linux-gnu-library/")
library(minet, lib.loc="~/R/x86_64-unknown-linux-gnu-library/")

# green tea
load("/protected/projects/pulmarray/Airway_Metanalysis/kyoto/aracne/data/grt_FS_MIM.rda")  # grt.mim
grt.clr.fs.network <- clr(grt.mim)
save(grt.clr.fs.network, file="/protected/projects/pulmarray/Airway_Metanalysis/kyoto/aracne/data/grt_FS_CLR.rda")
