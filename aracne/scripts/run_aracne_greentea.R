library(infotheo, lib.loc="~/R/x86_64-unknown-linux-gnu-library/")
library(minet, lib.loc="~/R/x86_64-unknown-linux-gnu-library/")

# green tea
lnames <- load("/protected/projects/pulmarray/Airway_Metanalysis/data/processed/input/clean/GreenTea_annotation_exprRMA_exprSCAN_clean_61.rda")
expr <- t(expr.rma)

grt.mim <- build.mim(expr, estimator="spearman")
save(grt.mim, file="/protected/projects/pulmarray/Airway_Metanalysis/kyoto/aracne/data/grt_MIM.rda")

grt.mim[1:5,1:5]
# load("/protected/projects/pulmarray/Airway_Metanalysis/kyoto/aracne/data/grt_MIM.rda")
grt.network <- aracne(grt.mim, eps=0)
save(grt.network, file="/protected/projects/pulmarray/Airway_Metanalysis/kyoto/aracne/data/grt_ARACNE.rda")
grt.network[1:5,1:5]