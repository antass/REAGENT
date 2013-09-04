library(infotheo, lib.loc="~/R/x86_64-unknown-linux-gnu-library/")
library(minet, lib.loc="~/R/x86_64-unknown-linux-gnu-library/")

# allegro
lnames <- load("/protected/projects/pulmarray/Airway_Metanalysis/data/processed/input/clean/Allegro16_annotation_expression_clean_ComBat_859.rda")
expr <- t(expression)

all.mim <- build.mim(expr, estimator="spearman")
save(all.mim, file="/protected/projects/pulmarray/Airway_Metanalysis/kyoto/aracne/data/allegro16_MIM.rda")
all.mim[1:5,1:5]

# load("/protected/projects/pulmarray/Airway_Metanalysis/kyoto/aracne/data/allegro_MIM.rda")
all.network <- aracne(all.mim, eps=0)
save(all.network, file="/protected/projects/pulmarray/Airway_Metanalysis/kyoto/aracne/data/allegro16_ARACNE.rda")
all.network[1:5,1:5]
