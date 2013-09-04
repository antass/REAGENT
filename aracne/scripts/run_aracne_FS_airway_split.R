library(infotheo, lib.loc="~/R/x86_64-unknown-linux-gnu-library/")
library(minet, lib.loc="~/R/x86_64-unknown-linux-gnu-library/")

# allegro
load("~/Meta_Analysis/kyoto/aracne/data/Feature_selection_DE_genes_cancer_premal_atleast2.rda")
lnames <- load("/protected/projects/pulmarray/Airway_Metanalysis/data/processed/input/clean/Airway_annotation_exprRMA_exprSCAN_clean_318.rda")
expr <- t(expr.rma)[, kyoto.de.genes.feat.sel.atleast2]

air.fs2.mim <- build.mim(expr, estimator="spearman")
save(air.fs2.mim, file="/protected/projects/pulmarray/Airway_Metanalysis/kyoto/aracne/data/airway_FS2_MIM.rda")
air.fs2.mim[1:5,1:5]

# load("/protected/projects/pulmarray/Airway_Metanalysis/kyoto/aracne/data/airegro_MIM.rda")
air.aracne.fs2.network <- aracne(air.fs2.mim, eps=0)
save(air.aracne.fs2.network, file="/protected/projects/pulmarray/Airway_Metanalysis/kyoto/aracne/data/airway_FS2_ARACNE.rda")
air.aracne.fs2.network[1:5,1:5]

air.clr.fs2.network <- clr(air.fs2.mim)
save(air.clr.fs2.network, file="/protected/projects/pulmarray/Airway_Metanalysis/kyoto/aracne/data/airway_FS2_CLR.rda")


# CANCER
ann.cancer <- ann[ !is.na(ann$cancer_status) & ann$cancer_status=="yes", ]
exp.cancer <- expr[rownames(expr) %in% ann.cancer$filename, ]

air.cancer.fs2.mim <- build.mim(exp.cancer, estimator="spearman")
save(air.cancer.fs2.mim, file="/protected/projects/pulmarray/Airway_Metanalysis/kyoto/aracne/data/airway_cancer_FS2_MIM.rda")
air.cancer.fs2.mim[1:5,1:5]

# load("/protected/projects/pulmarray/Airway_Metanalysis/kyoto/aracne/data/airegro_MIM.rda")
air.cancer.aracne.fs2.network <- aracne(air.cancer.fs2.mim, eps=0)
save(air.cancer.aracne.fs2.network, file="/protected/projects/pulmarray/Airway_Metanalysis/kyoto/aracne/data/airway_cancer_FS2_ARACNE.rda")
air.cancer.aracne.fs2.network[1:5,1:5]

air.cancer.clr.fs2.network <- clr(air.cancer.fs2.mim)
save(air.cancer.clr.fs2.network, file="/protected/projects/pulmarray/Airway_Metanalysis/kyoto/aracne/data/airway_cancer_FS2_CLR.rda")


# NORMAL
ann.normal <- ann[ !is.na(ann$cancer_status) & ann$cancer_status=="no", ]
exp.normal <- expr[rownames(expr) %in% ann.normal$filename, ]

air.normal.fs2.mim <- build.mim(exp.normal, estimator="spearman")
save(air.normal.fs2.mim, file="/protected/projects/pulmarray/Airway_Metanalysis/kyoto/aracne/data/airway_normal_FS2_MIM.rda")
air.normal.fs2.mim[1:5,1:5]

# load("/protected/projects/pulmarray/Airway_Metanalysis/kyoto/aracne/data/airegro_MIM.rda")
air.normal.aracne.fs2.network <- aracne(air.normal.fs2.mim, eps=0)
save(air.normal.aracne.fs2.network, file="/protected/projects/pulmarray/Airway_Metanalysis/kyoto/aracne/data/airway_normal_FS2_ARACNE.rda")
air.normal.aracne.fs2.network[1:5,1:5]

air.normal.clr.fs2.network <- clr(air.normal.fs2.mim)
save(air.normal.clr.fs2.network, file="/protected/projects/pulmarray/Airway_Metanalysis/kyoto/aracne/data/airway_normal_FS2_CLR.rda")
