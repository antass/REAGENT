library(infotheo, lib.loc="~/R/x86_64-unknown-linux-gnu-library/")
library(minet, lib.loc="~/R/x86_64-unknown-linux-gnu-library/")

load("~/Meta_Analysis/kyoto/aracne/data/Feature_selection_DE_genes_cancer_premal.rda")  # kyoto...
# allegro
lnames <- load("/protected/projects/pulmarray/Airway_Metanalysis/data/processed/input/clean/Allegro16_annotation_expression_clean_ComBat_859.rda")
expr <- t(expression)[, kyoto.de.genes.feat.sel]

all.fs.mim <- build.mim(expr, estimator="spearman")
save(all.fs.mim, file="/protected/projects/pulmarray/Airway_Metanalysis/kyoto/aracne/data/allegro16_FS_MIM.rda")
all.fs.mim[1:5,1:5]

# load("/protected/projects/pulmarray/Airway_Metanalysis/kyoto/aracne/data/allegro_MIM.rda")
all.aracne.fs.network <- aracne(all.fs.mim, eps=0)
save(all.aracne.fs.network, file="/protected/projects/pulmarray/Airway_Metanalysis/kyoto/aracne/data/allegro16_FS_ARACNE.rda")
all.aracne.fs.network[1:5,1:5]

all.clr.fs.network <- clr(all.fs.mim)
save(all.clr.fs.network, file="/protected/projects/pulmarray/Airway_Metanalysis/kyoto/aracne/data/allegro16_FS_CLR.rda")

