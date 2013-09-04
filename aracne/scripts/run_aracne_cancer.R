# library(infotheo, lib.loc="~/R/x86_64-unknown-linux-gnu-library/")
# library(minet, lib.loc="~/R/x86_64-unknown-linux-gnu-library/")
library(minet)

# joint
lnames <- load("~/Meta_Analysis/kyoto/aracne/data/Network_expression_cancer_936g_674s.rda")  # net.exp.c
expr <- t(net.exp.c)

joint.c.fs.mim <- build.mim(expr, estimator="spearman")
save(joint.c.fs.mim, file="~/Meta_Analysis/kyoto/aracne/data/joint_cancer_FS_MIM.rda")
# joint.c.fs.mim[1:5,1:5]

# load("/protected/projects/pulmarray/Airway_Metanalysis/kyoto/aracne/data/airegro_MIM.rda")
joint.aracne.c.fs.network <- aracne(joint.c.fs.mim, eps=0)
save(joint.aracne.c.fs.network, file="~/Meta_Analysis/kyoto/aracne/data/joint_cancer_FS_ARACNE.rda")
# joint.aracne.c.fs.network[1:5,1:5]

joint.clr.c.fs.network <- clr(joint.c.fs.mim)
save(joint.clr.c.fs.network, file="~/Meta_Analysis/kyoto/aracne/data/joint_cancer_FS_CLR.rda")
