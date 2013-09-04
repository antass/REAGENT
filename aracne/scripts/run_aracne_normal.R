# library(infotheo, lib.loc="~/R/x86_64-unknown-linux-gnu-library/")
# library(minet, lib.loc="~/R/x86_64-unknown-linux-gnu-library/")
library(minet)

# joint
lnames <- load("~/Meta_Analysis/kyoto/aracne/data/Network_expression_normal_936g_609s.rda")  # net.exp.n
expr <- t(net.exp.n)

joint.n.fs.mim <- build.mim(expr, estimator="spearman")
save(joint.n.fs.mim, file="~/Meta_Analysis/kyoto/aracne/data/joint_normal_FS_MIM.rda")
# joint.n.fs.mim[1:5,1:5]

# load("/protected/projects/pulmarray/Airway_Metanalysis/kyoto/aracne/data/airegro_MIM.rda")
joint.aracne.n.fs.network <- aracne(joint.n.fs.mim, eps=0)
save(joint.aracne.n.fs.network, file="~/Meta_Analysis/kyoto/aracne/data/joint_normal_FS_ARACNE.rda")
# joint.aracne.n.fs.network[1:5,1:5]

joint.clr.n.fs.network <- clr(joint.n.fs.mim)
save(joint.clr.n.fs.network, file="~/Meta_Analysis/kyoto/aracne/data/joint_normal_FS_CLR.rda")
