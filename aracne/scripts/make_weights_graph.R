### ARACNE
lname <- load("~/Meta_Analysis/kyoto/aracne/data/joint_normal_FS_ARACNE.rda")  # joint.aracne.n.fs.network
lname <- load("~/Meta_Analysis/kyoto/aracne/data/joint_cancer_FS_ARACNE.rda")  # joint.aracne.c.fs.network
lname <- load("~/Meta_Analysis/kyoto/aracne/data/joint_premalignancy_FS_ARACNE.rda")  # joint.aracne.p.fs.network
lname <- load("~/Meta_Analysis/kyoto/aracne/data/allegro_normal_FS_ARACNE.rda")  # all.aracne.n.fs.network
lname <- load("~/Meta_Analysis/kyoto/aracne/data/airway_normal_FS_ARACNE.rda")  # joint.aracne.c.fs.network
lname <- load("~/Meta_Analysis/kyoto/aracne/data/greentea_normal_FS_ARACNE.rda")  # joint.aracne.p.fs.networklname <- load("~/Meta_Analysis/kyoto/aracne/data/joint_normal_FS_ARACNE.rda")  # joint.aracne.n.fs.network
lname <- load("~/Meta_Analysis/kyoto/aracne/data/lam_normal_FS_ARACNE.rda")  # joint.aracne.c.fs.network
lname <- load("~/Meta_Analysis/kyoto/aracne/data/allegro_cancer_FS_ARACNE.rda")  # joint.aracne.p.fs.networklname <- load("~/Meta_Analysis/kyoto/aracne/data/joint_normal_FS_ARACNE.rda")  # joint.aracne.n.fs.network
lname <- load("~/Meta_Analysis/kyoto/aracne/data/airway_cancer_FS_ARACNE.rda")  # joint.aracne.c.fs.network
lname <- load("~/Meta_Analysis/kyoto/aracne/data/greentea_premalignancy_FS_ARACNE.rda")  # joint.aracne.p.fs.network
lname <- load("~/Meta_Analysis/kyoto/aracne/data/lam_premalignancy_FS_ARACNE.rda")  # joint.aracne.p.fs.network
### CLR
lname <- load("~/Meta_Analysis/kyoto/aracne/data/joint_normal_FS_CLR.rda")  # joint.aracne.n.fs.network
lname <- load("~/Meta_Analysis/kyoto/aracne/data/joint_cancer_FS_CLR.rda")  # joint.aracne.c.fs.network
lname <- load("~/Meta_Analysis/kyoto/aracne/data/joint_premalignancy_FS_CLR.rda")  # joint.aracne.p.fs.network
lname <- load("~/Meta_Analysis/kyoto/aracne/data/allegro_normal_FS_CLR.rda")  # all.aracne.n.fs.network
lname <- load("~/Meta_Analysis/kyoto/aracne/data/airway_normal_FS_CLR.rda")  # joint.aracne.c.fs.network
lname <- load("~/Meta_Analysis/kyoto/aracne/data/greentea_normal_FS_CLR.rda")  # joint.aracne.p.fs.networklname <- load("~/Meta_Analysis/kyoto/aracne/data/joint_normal_FS_CLR.rda")  # joint.aracne.n.fs.network
lname <- load("~/Meta_Analysis/kyoto/aracne/data/lam_normal_FS_CLR.rda")  # joint.aracne.c.fs.network
lname <- load("~/Meta_Analysis/kyoto/aracne/data/allegro_cancer_FS_CLR.rda")  # joint.aracne.p.fs.networklname <- load("~/Meta_Analysis/kyoto/aracne/data/joint_normal_FS_CLR.rda")  # joint.aracne.n.fs.network
lname <- load("~/Meta_Analysis/kyoto/aracne/data/airway_cancer_FS_CLR.rda")  # joint.aracne.c.fs.network
lname <- load("~/Meta_Analysis/kyoto/aracne/data/greentea_premalignancy_FS_CLR.rda")  # joint.aracne.p.fs.network
lname <- load("~/Meta_Analysis/kyoto/aracne/data/lam_premalignancy_FS_CLR.rda")  # joint.aracne.p.fs.network

n <- as.vector(joint.aracne.n.fs.network)
n.o <- n[order(n, decreasing=TRUE)]
n.o <- n.o[n.o>0]

c <- as.vector(joint.aracne.c.fs.network)
c.o <- c[order(c, decreasing=TRUE)]
c.o <- c.o[c.o>0]

p <- as.vector(joint.aracne.p.fs.network)
p.o <- p[order(p, decreasing=TRUE)]
p.o <- p.o[p.o>0]

xlim <- max(length(n.o), length(c.o), length(p.o))

plot(y=n.o, x=1:length(n.o), xlim=c(0,xlim), ylab="Weight", xlab="Edge weight order", type="l", col="green")
lines(y=c.o, x=1:length(c.o), xlim=c(0,xlim), col="red")
lines(y=p.o, x=1:length(p.o), xlim=c(0,xlim), col="blue")



n <- as.vector(joint.clr.n.fs.network)
n.o <- n[order(n, decreasing=TRUE)]
n.o <- n.o[n.o>0]

c <- as.vector(joint.clr.c.fs.network)
c.o <- c[order(c, decreasing=TRUE)]
c.o <- c.o[c.o>0]

p <- as.vector(joint.clr.p.fs.network)
p.o <- p[order(p, decreasing=TRUE)]
p.o <- p.o[p.o>0]

xlim <- max(length(n.o), length(c.o), length(p.o))

plot(y=n.o, x=1:length(n.o), xlim=c(0,xlim), ylab="Weight", xlab="Edge weight order", type="l", col="green")
lines(y=c.o, x=1:length(c.o), xlim=c(0,xlim), col="red")
lines(y=p.o, x=1:length(p.o), xlim=c(0,xlim), col="blue")

plot(y=n.o, x=1:length(n.o), xlim=c(0,1000), ylim=c(0,100), ylab="Weight", xlab="Edge weight order", type="l", col="green")
lines(y=c.o, x=1:length(c.o), xlim=c(0,xlim), col="red")
lines(y=p.o, x=1:length(p.o), xlim=c(0,xlim), col="blue")


