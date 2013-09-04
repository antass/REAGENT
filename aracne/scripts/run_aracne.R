library(infotheo, lib.loc="~/R/x86_64-unknown-linux-gnu-library/")

library(minet, lib.loc="~/R/x86_64-unknown-linux-gnu-library/")
# library(org.Hs.eg.db)

# # convert easily from affy probe ID to EntrezGene ID
# x <- org.Hs.egSYMBOL
# # Get the gene symbol that are mapped to an entrez gene identifiers
# mapped_genes <- mappedkeys(x)
# # Convert to a list
# xx <- as.list(x[mapped_genes])

# allegro
# lnames <- load("/protected/projects/pulmarray/Airway_Metanalysis/data/processed/input/clean/Allegro_annotation_exprRMA_exprSCAN_clean_ComBat_507.rda")
# expr <- t(expr.rma.combat)

# mim <- build.mim(expr, estimator="spearman")
# all.mim <- mim
# save(all.mim, file="/protected/projects/pulmarray/Airway_Metanalysis/kyoto/aracne/data/allegro_MIM.rda")

load("/protected/projects/pulmarray/Airway_Metanalysis/kyoto/aracne/data/allegro_MIM.rda")
mim <- all.mim
network <- aracne(mim, eps=0)
all.network <- network
save(all.network, file="/protected/projects/pulmarray/Airway_Metanalysis/kyoto/aracne/data/allegro_ARACNE.rda")


# airway
# lnames <- load("/protected/projects/pulmarray/Airway_Metanalysis/data/processed/input/clean/Airway_annotation_exprRMA_exprSCAN_clean_318.rda")
# expr <- t(expr.rma)
# 
# mim <- build.mim(expr, estimator="spearman")
# air.mim <- mim
# save(air.mim, file="/protected/projects/pulmarray/Airway_Metanalysis/kyoto/aracne/data/airway_MIM.rda")

load("/protected/projects/pulmarray/Airway_Metanalysis/kyoto/aracne/data/airway_MIM.rda")
mim <- air.mim
network <- aracne(mim, eps=0)
air.network <- network
save(air.network, file="/protected/projects/pulmarray/Airway_Metanalysis/kyoto/aracne/data/airway_ARACNE.rda")


# lam
# lnames <- load("/protected/projects/pulmarray/Airway_Metanalysis/data/processed/input/clean/Lam238_annotation_exprRMA_exprSCAN_clean_226.rda")
# expr <- t(expr.rma)
# 
# mim <- build.mim(expr, estimator="spearman")
# lam.mim <- mim
# save(lam.mim, file="/protected/projects/pulmarray/Airway_Metanalysis/kyoto/aracne/data/lam_MIM.rda")

load("/protected/projects/pulmarray/Airway_Metanalysis/kyoto/aracne/data/lam_MIM.rda")
mim <- lam.mim
network <- aracne(mim, eps=0)
lam.network <- network
save(lam.network, file="/protected/projects/pulmarray/Airway_Metanalysis/kyoto/aracne/data/lam_ARACNE.rda")


# green tea
lnames <- load("/protected/projects/pulmarray/Airway_Metanalysis/data/processed/input/clean/GreenTea_annotation_exprRMA_exprSCAN_clean_61.rda")
expr <- t(expr.rma)

mim <- build.mim(expr, estimator="spearman")
grt.mim <- mim
save(grt.mim, file="/protected/projects/pulmarray/Airway_Metanalysis/kyoto/aracne/data/grt_MIM.rda")

# load("/protected/projects/pulmarray/Airway_Metanalysis/kyoto/aracne/data/grt_MIM.rda")
# mim <- grt.mim
network <- aracne(mim, eps=0)
grt.network <- network
save(grt.network, file="/protected/projects/pulmarray/Airway_Metanalysis/kyoto/aracne/data/grt_ARACNE.rda")
