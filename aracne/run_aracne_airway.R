library(infotheo, lib.loc="~/R/x86_64-unknown-linux-gnu-library/")
library(minet, lib.loc="~/R/x86_64-unknown-linux-gnu-library/")

# airway
lnames <- load("/protected/projects/pulmarray/Airway_Metanalysis/data/processed/input/clean/Airway_annotation_exprRMA_exprSCAN_clean_318.rda")
expr <- t(expr.rma)

mim <- build.mim(expr, estimator="spearman")
air.mim <- mim
save(air.mim, file="/protected/projects/pulmarray/Airway_Metanalysis/kyoto/aracne/data/airway_MIM.rda")

load("/protected/projects/pulmarray/Airway_Metanalysis/kyoto/aracne/data/airway_MIM.rda")
mim <- air.mim
network <- aracne(mim, eps=0)
air.network <- network
save(air.network, file="/protected/projects/pulmarray/Airway_Metanalysis/kyoto/aracne/data/airway_ARACNE.rda")
