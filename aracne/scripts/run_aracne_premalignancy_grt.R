# library(infotheo, lib.loc="~/R/x86_64-unknown-linux-gnu-library/")
# library(minet, lib.loc="~/R/x86_64-unknown-linux-gnu-library/")
library(minet)


# greentea
lnames <- load("~/Meta_Analysis/kyoto/aracne/data/Network_expression_premalignancy_936g_120s.rda")  # grt.exp.ar.p
expr <- t(grt.exp.ar.p)

grt.p.fs.mim <- build.mim(expr, estimator="spearman")
save(grt.p.fs.mim, file="~/Meta_Analysis/kyoto/aracne/data/greentea_premalignancy_FS_MIM.rda")
# grt.p.fs.mim[1:5,1:5]

# load("/protected/projects/pulmarray/Airway_Metanalysis/kyoto/aracne/data/airegro_MIM.rda")
grt.aracne.p.fs.network <- aracne(grt.p.fs.mim, eps=0)
save(grt.aracne.p.fs.network, file="~/Meta_Analysis/kyoto/aracne/data/greentea_premalignancy_FS_ARACNE.rda")
# grt.aracne.p.fs.network[1:5,1:5]

grt.clr.p.fs.network <- clr(grt.p.fs.mim)
save(grt.clr.p.fs.network, file="~/Meta_Analysis/kyoto/aracne/data/greentea_premalignancy_FS_CLR.rda")
