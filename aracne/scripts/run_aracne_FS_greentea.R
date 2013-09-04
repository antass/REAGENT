library(infotheo, lib.loc="~/R/x86_64-unknown-linux-gnu-library/")
library(minet, lib.loc="~/R/x86_64-unknown-linux-gnu-library/")

# allegro
load("~/Meta_Analysis/kyoto/aracne/data/Feature_selection_DE_genes_cancer_premal.rda")  # kyoto...
lnames <- load("/protected/projects/pulmarray/Airway_Metanalysis/data/processed/input/clean/GreenTea_annotation_premalignancy_exprRMA_exprSCAN_clean_61.rda")
expr <- t(expr.rma)[, kyoto.de.genes.feat.sel]

grt.mim <- build.mim(expr, estimator="spearman")
save(grt.mim, file="/protected/projects/pulmarray/Airway_Metanalysis/kyoto/aracne/data/grt_FS_MIM.rda")
grt.mim[1:5,1:5]

# load("/protected/projects/pulmarray/grtway_Metanalysis/kyoto/aracne/data/grtegro_MIM.rda")
grt.aracne.fs.network <- aracne(grt.mim, eps=0)
save(grt.aracne.fs.network, file="/protected/projects/pulmarray/Airway_Metanalysis/kyoto/aracne/data/grt_FS_ARACNE.rda")
grt.aracne.fs.network[1:5,1:5]

grt.clr.fs.network <- clr(grt.mim)
save(grt.clr.fs.network, file="/protected/projects/pulmarray/Airway_Metanalysis/kyoto/aracne/data/grt_FS_CLR.rda")

