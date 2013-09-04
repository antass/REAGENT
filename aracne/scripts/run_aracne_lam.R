library(infotheo, lib.loc="~/R/x86_64-unknown-linux-gnu-library/")
library(minet, lib.loc="~/R/x86_64-unknown-linux-gnu-library/")

# lam
lnames <- load("/protected/projects/pulmarray/Airway_Metanalysis/data/processed/input/clean/Lam238_annotation_exprRMA_exprSCAN_clean_226.rda")
expr <- t(expr.rma)

lam.mim <- build.mim(expr, estimator="spearman")
save(lam.mim, file="/protected/projects/pulmarray/Airway_Metanalysis/kyoto/aracne/data/lam_MIM.rda")


lam.mim[1:5,1:5]

load("/protected/projects/pulmarray/Airway_Metanalysis/kyoto/aracne/data/lam_MIM.rda")
lam.network <- aracne(lam.mim, eps=0)
save(lam.network, file="/protected/projects/pulmarray/Airway_Metanalysis/kyoto/aracne/data/lam_ARACNE.rda")
lam.network[1:5,1:5]