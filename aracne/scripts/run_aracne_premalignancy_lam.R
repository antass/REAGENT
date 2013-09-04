# library(infotheo, lib.loc="~/R/x86_64-unknown-linux-gnu-library/")
# library(minet, lib.loc="~/R/x86_64-unknown-linux-gnu-library/")
library(minet)


# lam
lnames <- load("~/Meta_Analysis/kyoto/aracne/data/Network_expression_premalignancy_lam_936g_71s.rda")  # lam.exp.ar.p
expr <- t(lam.exp.ar.p)

lam.p.fs.mim <- build.mim(expr, estimator="spearman")
save(lam.p.fs.mim, file="~/Meta_Analysis/kyoto/aracne/data/lam_premalignancy_FS_MIM.rda")
# lam.p.fs.mim[1:5,1:5]

# load("/protected/projects/pulmarray/Airway_Metanalysis/kyoto/aracne/data/airegro_MIM.rda")
lam.aracne.p.fs.network <- aracne(lam.p.fs.mim, eps=0)
save(lam.aracne.p.fs.network, file="~/Meta_Analysis/kyoto/aracne/data/lam_premalignancy_FS_ARACNE.rda")
# lam.aracne.p.fs.network[1:5,1:5]

lam.clr.p.fs.network <- clr(lam.p.fs.mim)
save(lam.clr.p.fs.network, file="~/Meta_Analysis/kyoto/aracne/data/lam_premalignancy_FS_CLR.rda")
