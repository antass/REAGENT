# library(infotheo, lib.loc="~/R/x86_64-unknown-linux-gnu-library/")
# library(minet, lib.loc="~/R/x86_64-unknown-linux-gnu-library/")
library(minet)

# allegro
lnames <- load("~/Meta_Analysis/kyoto/aracne/data/Network_expression_cancer_allegro_936g_555s.rda")  # all.exp.ar.c
expr <- t(all.exp.ar.c)

all.c.fs.mim <- build.mim(expr, estimator="spearman")
save(all.c.fs.mim, file="~/Meta_Analysis/kyoto/aracne/data/allegro_cancer_FS_MIM.rda")
# all.c.fs.mim[1:5,1:5]

# load("/protected/projects/pulmarray/Airway_Metanalysis/kyoto/aracne/data/airegro_MIM.rda")
all.aracne.c.fs.network <- aracne(all.c.fs.mim, eps=0)
save(all.aracne.c.fs.network, file="~/Meta_Analysis/kyoto/aracne/data/allegro_cancer_FS_ARACNE.rda")
# all.aracne.c.fs.network[1:5,1:5]

all.clr.c.fs.network <- clr(all.c.fs.mim)
save(all.clr.c.fs.network, file="~/Meta_Analysis/kyoto/aracne/data/allegro_cancer_FS_CLR.rda")
