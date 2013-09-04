library(infotheo, lib.loc="~/R/x86_64-unknown-linux-gnu-library/")
library(minet, lib.loc="~/R/x86_64-unknown-linux-gnu-library/")

load("~/Meta_Analysis/kyoto/aracne/data/Feature_selection_DE_genes_cancer_premal_atleast2.rda")  # kyoto...
# grt
lnames <- load("/protected/projects/pulmarray/Airway_Metanalysis/data/processed/input/clean/GreenTea_annotation_premalignancy_exprRMA_exprSCAN_clean_61.rda")
expr <- t(expr.rma)[, kyoto.de.genes.feat.sel.atleast2]

grt.fs2.mim <- build.mim(expr, estimator="spearman")
save(grt.fs2.mim, file="/protected/projects/pulmarray/Airway_Metanalysis/kyoto/aracne/data/grt_FS2_MIM.rda")
grt.fs2.mim[1:5,1:5]

# load("/protected/projects/pulmarray/Airway_Metanalysis/kyoto/aracne/data/airegro_MIM.rda")
grt.aracne.fs2.network <- aracne(grt.fs2.mim, eps=0)
save(grt.aracne.fs2.network, file="/protected/projects/pulmarray/Airway_Metanalysis/kyoto/aracne/data/grt_FS2_ARACNE.rda")
grt.aracne.fs2.network[1:5,1:5]

grt.clr.fs2.network <- clr(grt.fs2.mim)
save(grt.clr.fs2.network, file="/protected/projects/pulmarray/Airway_Metanalysis/kyoto/aracne/data/grt_FS2_CLR.rda")

# PREMALIGNANT LESIONS
ann.dyspl <- ann.premal[ !is.na(ann.premal$DYSGR_CAT) & !(ann.premal$DYSGR_CAT%in%c("Normal", "CarcinomaOther", "NoData", "Unsatis")), ]
exp.dyspl <- expr[rownames(expr) %in% ann.dyspl$filename, ]

grt.premal.fs2.mim <- build.mim(exp.dyspl, estimator="spearman")
save(grt.premal.fs2.mim, file="/protected/projects/pulmarray/Airway_Metanalysis/kyoto/aracne/data/grt_premal_FS2_MIM.rda")
grt.premal.fs2.mim[1:5,1:5]

# load("/protected/projects/pulmarray/Airway_Metanalysis/kyoto/aracne/data/airegro_MIM.rda")
grt.premal.aracne.fs2.network <- aracne(grt.premal.fs2.mim, eps=0)
save(grt.premal.aracne.fs2.network, file="/protected/projects/pulmarray/Airway_Metanalysis/kyoto/aracne/data/grt_premal_FS2_ARACNE.rda")
grt.premal.aracne.fs2.network[1:5,1:5]
dim(grt.premal.aracne.fs2.network)

grt.premal.clr.fs2.network <- clr(grt.premal.fs2.mim)
save(grt.premal.clr.fs2.network, file="/protected/projects/pulmarray/Airway_Metanalysis/kyoto/aracne/data/grt_premal_FS2_CLR.rda")
dim(grt.premal.aracne.fs2.network)


# NORMAL
ann.normal <- ann.premal[ !is.na(ann.premal$DYSGR_CAT) & !(ann.premal$DYSGR_CAT%in%c("Normal", "CarcinomaOther", "NoData", "Unsatis")), ]
exp.normal <- expr[rownames(expr) %in% ann.normal$filename, ]

grt.normal.fs2.mim <- build.mim(exp.normal, estimator="spearman")
save(grt.normal.fs2.mim, file="/protected/projects/pulmarray/Airway_Metanalysis/kyoto/aracne/data/grt_normal_FS2_MIM.rda")
grt.normal.fs2.mim[1:5,1:5]

# load("/protected/projects/pulmarray/Airway_Metanalysis/kyoto/aracne/data/airegro_MIM.rda")
grt.normal.aracne.fs2.network <- aracne(grt.normal.fs2.mim, eps=0)
save(grt.normal.aracne.fs2.network, file="/protected/projects/pulmarray/Airway_Metanalysis/kyoto/aracne/data/grt_normal_FS2_ARACNE.rda")
grt.normal.aracne.fs2.network[1:5,1:5]

grt.normal.clr.fs2.network <- clr(grt.normal.fs2.mim)
save(grt.normal.clr.fs2.network, file="/protected/projects/pulmarray/Airway_Metanalysis/kyoto/aracne/data/grt_normal_FS2_CLR.rda")


