# library(infotheo, lib.loc="~/R/x86_64-unknown-linux-gnu-library/")
# library(minet, lib.loc="~/R/x86_64-unknown-linux-gnu-library/")
library(minet)


# joint
lnames <- load("~/Meta_Analysis/kyoto/aracne/data/Network_expression_premalignancy_936g_120s.rda")  # net.exp.p
expr <- t(net.exp.p)

joint.p.fs.mim <- build.mim(expr, estimator="spearman")
save(joint.p.fs.mim, file="~/Meta_Analysis/kyoto/aracne/data/joint_premalignancy_FS_MIM.rda")
# joint.p.fs.mim[1:5,1:5]

# load("/protected/projects/pulmarray/Airway_Metanalysis/kyoto/aracne/data/airegro_MIM.rda")
joint.aracne.p.fs.network <- aracne(joint.p.fs.mim, eps=0)
save(joint.aracne.p.fs.network, file="~/Meta_Analysis/kyoto/aracne/data/joint_premalignancy_FS_ARACNE.rda")
# joint.aracne.p.fs.network[1:5,1:5]

joint.clr.p.fs.network <- clr(joint.p.fs.mim)
save(joint.clr.p.fs.network, file="~/Meta_Analysis/kyoto/aracne/data/joint_premalignancy_FS_CLR.rda")
