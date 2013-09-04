# library(infotheo, lib.loc="~/R/x86_64-unknown-linux-gnu-library/")
# library(minet, lib.loc="~/R/x86_64-unknown-linux-gnu-library/")
library(minet)

# greentea
lnames <- load("~/Meta_Analysis/kyoto/aracne/data/Network_expression_normal_greentea_936g_11s.rda")  # grt.exp.ar.n
expr <- t(grt.exp.ar.n)

grt.n.fs.mim <- build.mim(expr, estimator="spearman")
save(grt.n.fs.mim, file="~/Meta_Analysis/kyoto/aracne/data/greentea_normal_FS_MIM.rda")
# grt.n.fs.mim[1:5,1:5]

# load("/protected/projects/pulmarray/Airway_Metanalysis/kyoto/aracne/data/airegro_MIM.rda")
grt.aracne.n.fs.network <- aracne(grt.n.fs.mim, eps=0)
save(grt.aracne.n.fs.network, file="~/Meta_Analysis/kyoto/aracne/data/greentea_normal_FS_ARACNE.rda")
# grt.aracne.n.fs.network[1:5,1:5]

grt.clr.n.fs.network <- clr(grt.n.fs.mim)
save(grt.clr.n.fs.network, file="~/Meta_Analysis/kyoto/aracne/data/greentea_normal_FS_CLR.rda")
