library(infotheo, lib.loc="~/R/x86_64-unknown-linux-gnu-library/")
library(minet, lib.loc="~/R/x86_64-unknown-linux-gnu-library/")

# allegro
load("~/Meta_Analysis/kyoto/aracne/data/Feature_selection_DE_genes_cancer_premal.rda")
lnames <- load("/protected/projects/pulmarray/Airway_Metanalysis/data/processed/input/clean/Airway_annotation_exprRMA_exprSCAN_clean_318.rda")
expr <- t(expr.rma)[, kyoto.de.genes.feat.sel]

air.fs.mim <- build.mim(expr, estimator="spearman")
save(air.fs.mim, file="/protected/projects/pulmarray/Airway_Metanalysis/kyoto/aracne/data/airway_FS_MIM.rda")
air.fs.mim[1:5,1:5]

# load("/protected/projects/pulmarray/Airway_Metanalysis/kyoto/aracne/data/airegro_MIM.rda")
air.aracne.fs.network <- aracne(air.fs.mim, eps=0)
save(air.aracne.fs.network, file="/protected/projects/pulmarray/Airway_Metanalysis/kyoto/aracne/data/airway_FS_ARACNE.rda")
air.aracne.fs.network[1:5,1:5]

air.clr.fs.network <- clr(air.fs.mim)
save(air.clr.fs.network, file="/protected/projects/pulmarray/Airway_Metanalysis/kyoto/aracne/data/airway_FS_CLR.rda")
