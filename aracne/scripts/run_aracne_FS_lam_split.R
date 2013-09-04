library(infotheo, lib.loc="~/R/x86_64-unknown-linux-gnu-library/")
library(minet, lib.loc="~/R/x86_64-unknown-linux-gnu-library/")

load("~/Meta_Analysis/kyoto/aracne/data/Feature_selection_DE_genes_cancer_premal_atleast2.rda")  # kyoto...
# lam
lnames <- load("/protected/projects/pulmarray/Airway_Metanalysis/data/processed/input/clean/Lam238_annotation_premalignancy_exprRMA_exprSCAN_clean_226.rda")
expr <- t(expr.rma)[, kyoto.de.genes.feat.sel.atleast2]

lam.fs2.mim <- build.mim(expr, estimator="spearman")
save(lam.fs2.mim, file="/protected/projects/pulmarray/Airway_Metanalysis/kyoto/aracne/data/lam_FS2_MIM.rda")
lam.fs2.mim[1:5,1:5]

# load("/protected/projects/pulmarray/Airway_Metanalysis/kyoto/aracne/data/airegro_MIM.rda")
lam.aracne.fs2.network <- aracne(lam.fs2.mim, eps=0)
save(lam.aracne.fs2.network, file="/protected/projects/pulmarray/Airway_Metanalysis/kyoto/aracne/data/lam_FS2_ARACNE.rda")
lam.aracne.fs2.network[1:5,1:5]

lam.clr.fs2.network <- clr(lam.fs2.mim)
save(lam.clr.fs2.network, file="/protected/projects/pulmarray/Airway_Metanalysis/kyoto/aracne/data/lam_FS2_CLR.rda")

# PREMALIGNANT LESIONS
ann.dyspl <- ann.premal[ !is.na(ann.premal$DYSGR_CAT) & !(ann.premal$DYSGR_CAT%in%c("Normal", "CarcinomaOther", "NoData", "Unsatis")), ]
exp.dyspl <- expr[rownames(expr) %in% ann.dyspl$filename, ]

lam.premal.fs2.mim <- build.mim(exp.dyspl, estimator="spearman")
save(lam.premal.fs2.mim, file="/protected/projects/pulmarray/Airway_Metanalysis/kyoto/aracne/data/lam_premal_FS2_MIM.rda")
lam.premal.fs2.mim[1:5,1:5]

# load("/protected/projects/pulmarray/Airway_Metanalysis/kyoto/aracne/data/airegro_MIM.rda")
lam.premal.aracne.fs2.network <- aracne(lam.premal.fs2.mim, eps=0)
save(lam.premal.aracne.fs2.network, file="/protected/projects/pulmarray/Airway_Metanalysis/kyoto/aracne/data/lam_premal_FS2_ARACNE.rda")
lam.premal.aracne.fs2.network[1:5,1:5]
dim(lam.premal.aracne.fs2.network)

lam.premal.clr.fs2.network <- clr(lam.premal.fs2.mim)
save(lam.premal.clr.fs2.network, file="/protected/projects/pulmarray/Airway_Metanalysis/kyoto/aracne/data/lam_premal_FS2_CLR.rda")
dim(lam.premal.aracne.fs2.network)


# NORMAL
ann.normal <- ann.premal[ !is.na(ann.premal$DYSGR_CAT) & !(ann.premal$DYSGR_CAT%in%c("Normal", "CarcinomaOther", "NoData", "Unsatis")), ]
exp.normal <- expr[rownames(expr) %in% ann.normal$filename, ]

lam.normal.fs2.mim <- build.mim(exp.normal, estimator="spearman")
save(lam.normal.fs2.mim, file="/protected/projects/pulmarray/Airway_Metanalysis/kyoto/aracne/data/lam_normal_FS2_MIM.rda")
lam.normal.fs2.mim[1:5,1:5]

# load("/protected/projects/pulmarray/Airway_Metanalysis/kyoto/aracne/data/airegro_MIM.rda")
lam.normal.aracne.fs2.network <- aracne(lam.normal.fs2.mim, eps=0)
save(lam.normal.aracne.fs2.network, file="/protected/projects/pulmarray/Airway_Metanalysis/kyoto/aracne/data/lam_normal_FS2_ARACNE.rda")
lam.normal.aracne.fs2.network[1:5,1:5]

lam.normal.clr.fs2.network <- clr(lam.normal.fs2.mim)
save(lam.normal.clr.fs2.network, file="/protected/projects/pulmarray/Airway_Metanalysis/kyoto/aracne/data/lam_normal_FS2_CLR.rda")


