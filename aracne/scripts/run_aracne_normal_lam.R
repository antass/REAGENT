# library(infotheo, lib.loc="~/R/x86_64-unknown-linux-gnu-library/")
# library(minet, lib.loc="~/R/x86_64-unknown-linux-gnu-library/")
library(minet)

# joint
lnames <- load("~/Meta_Analysis/kyoto/aracne/data/Network_expression_normal_lam_936g_147s.rda")  # net.exp.n
expr <- t(lam.exp.ar.n)

lam.n.fs.mim <- build.mim(expr, estimator="spearman")
save(lam.n.fs.mim, file="~/Meta_Analysis/kyoto/aracne/data/lam_normal_FS_MIM.rda")
# lam.n.fs.mim[1:5,1:5]

# load("/protected/projects/pulmarray/Airway_Metanalysis/kyoto/aracne/data/airegro_MIM.rda")
lam.aracne.n.fs.network <- aracne(lam.n.fs.mim, eps=0)
save(lam.aracne.n.fs.network, file="~/Meta_Analysis/kyoto/aracne/data/lam_normal_FS_ARACNE.rda")
# lam.aracne.n.fs.network[1:5,1:5]

lam.clr.n.fs.network <- clr(lam.n.fs.mim)
save(lam.clr.n.fs.network, file="~/Meta_Analysis/kyoto/aracne/data/lam_normal_FS_CLR.rda")
