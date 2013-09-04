# Using 2of3 strategy for genes
genes <- kyoto.de.genes.feat.sel.atleast2

grt.exp.ar <- grt.exp[genes, ]
all.exp.ar <- all.exp[genes, ]
air.exp.ar <- air.exp[genes, ]
lam.exp.ar <- lam.exp[genes, ]

grt.ann.ar <- grt.ann
grt.ann.ar$DYSGR_CAT <- as.character(grt.ann.ar$DYSGR_CAT)
grt.ann.ar$DYSGR_CAT[ grt.ann.ar$DYSGR_CAT %in% c("CarcinomaOther", "Unsatis")] <- "Other"
# grt.ann.ar$DYSGR_CAT[ grt.ann.ar$DYSGR_CAT == "NoData"] <- NA
grt.ann.ar$DYSGR_CAT[ grt.ann.ar$DYSGR_CAT == "Hyperplasia"] <- "Normal"
grt.ann.ar$DYSGR_CAT[ !grt.ann.ar$DYSGR_CAT %in% c("Normal", "Other")] <- "Dysplasia"
grt.ann.ar$DYSGR_CAT <- droplevels(as.factor(grt.ann.ar$DYSGR_CAT))
grt.ann$DYSGR_CAT3 <- grt.ann.ar$DYSGR_CAT

lam.ann.ar <- lam.ann
lam.ann.ar$DYSGR_CAT <- as.character(lam.ann.ar$DYSGR_CAT)
lam.ann.ar$DYSGR_CAT[ lam.ann.ar$DYSGR_CAT %in% c("CarcinomaOther", "Unsatis")] <- "Other"
lam.ann.ar$DYSGR_CAT[ lam.ann.ar$DYSGR_CAT == "NoData"] <- NA
lam.ann.ar$DYSGR_CAT[ lam.ann.ar$DYSGR_CAT == "Hyperplasia"] <- "Normal"
lam.ann.ar$DYSGR_CAT[ !is.na(lam.ann.ar$DYSGR_CAT) & !lam.ann.ar$DYSGR_CAT %in% c("Normal", "Other", NA)] <- "Dysplasia"
lam.ann.ar$DYSGR_CAT <- droplevels(as.factor(lam.ann.ar$DYSGR_CAT))
lam.ann$DYSGR_CAT3 <- lam.ann.ar$DYSGR_CAT


# NORMAL network (all, air, lam, grt)
grt.exp.ar.n <- grt.exp.ar[, which(grt.ann.ar$DYSGR_CAT=="Normal")]
all.exp.ar.n <- all.exp.ar[, which(all.ann$cancer_status=="no cancer")]
air.exp.ar.n <- air.exp.ar[, which(air.ann$cancer_status=="no")]
lam.exp.ar.n <- lam.exp.ar[, which(lam.ann.ar$DYSGR_CAT=="Normal")]
net.exp.n <- cbind(air.exp.ar.n, all.exp.ar.n, grt.exp.ar.n, lam.exp.ar.n)
save(net.exp.n, file="~/Meta_Analysis/kyoto/aracne/data/Network_expression_normal_936g_609s.rda")
save(air.exp.ar.n, file="~/Meta_Analysis/kyoto/aracne/data/Network_expression_normal_airway_936g_149s.rda")
save(all.exp.ar.n, file="~/Meta_Analysis/kyoto/aracne/data/Network_expression_normal_allegro_936g_302s.rda")
save(grt.exp.ar.n, file="~/Meta_Analysis/kyoto/aracne/data/Network_expression_normal_greentea_936g_11s.rda")
save(lam.exp.ar.n, file="~/Meta_Analysis/kyoto/aracne/data/Network_expression_normal_lam_936g_147s.rda")


# CANCER network (all, air)
all.exp.ar.c <- all.exp.ar[, which(all.ann$cancer_status=="cancer")]
air.exp.ar.c <- air.exp.ar[, which(air.ann$cancer_status=="yes")]
net.exp.c <- cbind(air.exp.ar.c, all.exp.ar.c)
save(net.exp.c, file="~/Meta_Analysis/kyoto/aracne/data/Network_expression_cancer_936g_674s.rda")
save(air.exp.ar.c, file="~/Meta_Analysis/kyoto/aracne/data/Network_expression_cancer_airway_936g_119s.rda")
save(all.exp.ar.c, file="~/Meta_Analysis/kyoto/aracne/data/Network_expression_cancer_allegro_936g_555s.rda")



# PREMALIGNANCY network (lam, grt)
grt.exp.ar.p <- grt.exp.ar[, which(grt.ann.ar$DYSGR_CAT=="Dysplasia")]
lam.exp.ar.p <- lam.exp.ar[, which(lam.ann.ar$DYSGR_CAT=="Dysplasia")]
net.exp.p <- cbind(grt.exp.ar.p, lam.exp.ar.p)
save(net.exp.p, file="~/Meta_Analysis/kyoto/aracne/data/Network_expression_premalignancy_936g_120s.rda")
save(grt.exp.ar.p, file="~/Meta_Analysis/kyoto/aracne/data/Network_expression_premalignancy_greentea_936g_49s.rda")
save(lam.exp.ar.p, file="~/Meta_Analysis/kyoto/aracne/data/Network_expression_premalignancy_lam_936g_71s.rda")





