# library(infotheo, lib.loc="~/R/x86_64-unknown-linux-gnu-library/")
# library(minet, lib.loc="~/R/x86_64-unknown-linux-gnu-library/")
library(minet)

# allegro
lnames <- load("~/Meta_Analysis/kyoto/aracne/data/Network_expression_normal_allegro_936g_302s.rda")  # all.exp.ar.n
expr <- t(all.exp.ar.n)

all.n.fs.mim <- build.mim(expr, estimator="spearman")
save(all.n.fs.mim, file="~/Meta_Analysis/kyoto/aracne/data/allegro_normal_FS_MIM.rda")
# all.n.fs.mim[1:5,1:5]

# load("/protected/projects/pulmarray/Airway_Metanalysis/kyoto/aracne/data/airegro_MIM.rda")
all.aracne.n.fs.network <- aracne(all.n.fs.mim, eps=0)
save(all.aracne.n.fs.network, file="~/Meta_Analysis/kyoto/aracne/data/allegro_normal_FS_ARACNE.rda")
# all.aracne.n.fs.network[1:5,1:5]

all.clr.n.fs.network <- clr(all.n.fs.mim)
save(all.clr.n.fs.network, file="~/Meta_Analysis/kyoto/aracne/data/allegro_normal_FS_CLR.rda")
