library(infotheo, lib.loc="~/R/x86_64-unknown-linux-gnu-library/")
library(minet, lib.loc="~/R/x86_64-unknown-linux-gnu-library/")

# allegro
lnames <- load("/protected/projects/pulmarray/Airway_Metanalysis/data/processed/input/clean/Allegro_annotation_exprRMA_exprSCAN_clean_ComBat_507.rda")
expr <- t(expr.rma.combat)

mim <- build.mim(expr, estimator="spearman")
all.mim <- mim
save(all.mim, file="/protected/projects/pulmarray/Airway_Metanalysis/kyoto/aracne/data/allegro_MIM.rda")

load("/protected/projects/pulmarray/Airway_Metanalysis/kyoto/aracne/data/allegro_MIM.rda")
mim <- all.mim
network <- aracne(mim, eps=0)
all.network <- network
save(all.network, file="/protected/projects/pulmarray/Airway_Metanalysis/kyoto/aracne/data/allegro_ARACNE.rda")

