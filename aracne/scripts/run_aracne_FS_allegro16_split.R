library(infotheo, lib.loc="~/R/x86_64-unknown-linux-gnu-library/")
library(minet, lib.loc="~/R/x86_64-unknown-linux-gnu-library/")

load("~/Meta_Analysis/kyoto/aracne/data/Feature_selection_DE_genes_cancer_premal_atleast2.rda")  # kyoto...
# allegro
lnames <- load("/protected/projects/pulmarray/Airway_Metanalysis/data/processed/input/clean/Allegro16_annotation_expression_clean_ComBat_859.rda")
expr <- t(expression)[, kyoto.de.genes.feat.sel.atleast2]
ann <- annotation

all.fs2.mim <- build.mim(expr, estimator="spearman")
save(all.fs2.mim, file="/protected/projects/pulmarray/Airway_Metanalysis/kyoto/aracne/data/allegro_FS2_MIM.rda")
all.fs2.mim[1:5,1:5]

# load("/protected/projects/pulmarray/Airway_Metanalysis/kyoto/aracne/data/airegro_MIM.rda")
all.aracne.fs2.network <- aracne(all.fs2.mim, eps=0)
save(all.aracne.fs2.network, file="/protected/projects/pulmarray/Airway_Metanalysis/kyoto/aracne/data/allegro_FS2_ARACNE.rda")
all.aracne.fs2.network[1:5,1:5]

all.clr.fs2.network <- clr(all.fs2.mim)
save(all.clr.fs2.network, file="/protected/projects/pulmarray/Airway_Metanalysis/kyoto/aracne/data/allegro_FS2_CLR.rda")

# CANCER
ann.cancer <- ann[ !is.na(ann$cancer_status) & ann$cancer_status=="yes", ]
exp.cancer <- expr[rownames(expr) %in% ann.cancer$filename, ]

all.cancer.fs2.mim <- build.mim(exp.cancer, estimator="spearman")
save(all.cancer.fs2.mim, file="/protected/projects/pulmarray/Airway_Metanalysis/kyoto/aracne/data/allegro_cancer_FS2_MIM.rda")
all.cancer.fs2.mim[1:5,1:5]

# load("/protected/projects/pulmarray/Airway_Metanalysis/kyoto/aracne/data/airegro_MIM.rda")
all.cancer.aracne.fs2.network <- aracne(all.cancer.fs2.mim, eps=0)
save(all.cancer.aracne.fs2.network, file="/protected/projects/pulmarray/Airway_Metanalysis/kyoto/aracne/data/allegro_cancer_FS2_ARACNE.rda")
all.cancer.aracne.fs2.network[1:5,1:5]

all.cancer.clr.fs2.network <- clr(all.cancer.fs2.mim)
save(all.cancer.clr.fs2.network, file="/protected/projects/pulmarray/Airway_Metanalysis/kyoto/aracne/data/allegro_cancer_FS2_CLR.rda")


# NORMAL
ann.normal <- ann[ !is.na(ann$cancer_status) & ann$cancer_status=="no", ]
exp.normal <- expr[rownames(expr) %in% ann.normal$filename, ]

all.normal.fs2.mim <- build.mim(exp.normal, estimator="spearman")
save(all.normal.fs2.mim, file="/protected/projects/pulmarray/Airway_Metanalysis/kyoto/aracne/data/allegro_normal_FS2_MIM.rda")
all.normal.fs2.mim[1:5,1:5]

# load("/protected/projects/pulmarray/Airway_Metanalysis/kyoto/aracne/data/airegro_MIM.rda")
all.normal.aracne.fs2.network <- aracne(all.normal.fs2.mim, eps=0)
save(all.normal.aracne.fs2.network, file="/protected/projects/pulmarray/Airway_Metanalysis/kyoto/aracne/data/allegro_normal_FS2_ARACNE.rda")
all.normal.aracne.fs2.network[1:5,1:5]

all.normal.clr.fs2.network <- clr(all.normal.fs2.mim)
save(all.normal.clr.fs2.network, file="/protected/projects/pulmarray/Airway_Metanalysis/kyoto/aracne/data/allegro_normal_FS2_CLR.rda")


