# library(infotheo, lib.loc="~/R/x86_64-unknown-linux-gnu-library/")
# library(minet, lib.loc="~/R/x86_64-unknown-linux-gnu-library/")
library(minet)

# airway
lnames <- load("~/Meta_Analysis/kyoto/aracne/data/Network_expression_cancer_airway_936g_119s.rda")  # air.exp.ar.c
expr <- t(air.exp.ar.c)

air.c.fs.mim <- build.mim(expr, estimator="spearman")
save(air.c.fs.mim, file="~/Meta_Analysis/kyoto/aracne/data/airway_cancer_FS_MIM.rda")
# air.c.fs.mim[1:5,1:5]

# load("/protected/projects/pulmarray/Airway_Metanalysis/kyoto/aracne/data/airegro_MIM.rda")
air.aracne.c.fs.network <- aracne(air.c.fs.mim, eps=0)
save(air.aracne.c.fs.network, file="~/Meta_Analysis/kyoto/aracne/data/airway_cancer_FS_ARACNE.rda")
# air.aracne.c.fs.network[1:5,1:5]

air.clr.c.fs.network <- clr(air.c.fs.mim)
save(air.clr.c.fs.network, file="~/Meta_Analysis/kyoto/aracne/data/airway_cancer_FS_CLR.rda")
