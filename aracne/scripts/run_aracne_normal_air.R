# library(infotheo, lib.loc="~/R/x86_64-unknown-linux-gnu-library/")
# library(minet, lib.loc="~/R/x86_64-unknown-linux-gnu-library/")
library(minet)

# airway
lnames <- load("~/Meta_Analysis/kyoto/aracne/data/Network_expression_normal_airway_936g_149s.rda")  # air.exp.ar.n
expr <- t(air.exp.ar.n)

air.n.fs.mim <- build.mim(expr, estimator="spearman")
save(air.n.fs.mim, file="~/Meta_Analysis/kyoto/aracne/data/airway_normal_FS_MIM.rda")
# air.n.fs.mim[1:5,1:5]

# load("/protected/projects/pulmarray/Airway_Metanalysis/kyoto/aracne/data/airegro_MIM.rda")
air.aracne.n.fs.network <- aracne(air.n.fs.mim, eps=0)
save(air.aracne.n.fs.network, file="~/Meta_Analysis/kyoto/aracne/data/airway_normal_FS_ARACNE.rda")
# air.aracne.n.fs.network[1:5,1:5]

air.clr.n.fs.network <- clr(air.n.fs.mim)
save(air.clr.n.fs.network, file="~/Meta_Analysis/kyoto/aracne/data/airway_normal_FS_CLR.rda")
