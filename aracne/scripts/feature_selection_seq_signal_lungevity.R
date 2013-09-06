# library(limma, lib.loc="~/R/x86_64-unknown-linux-gnu-library/")
library(limma)

load("~/Meta_Analysis/kyoto/aracne/data/exp_common_genes_all_air_lam_grt_lgv_8239g.rda")  # all.exp.common ...

annotation <- lgv.ann
expression <- lgv.exp.common 

annot.base <- annotation

data.0 <- expression[, as.character(annot.base$filename)]
data.0 <- data.0[, match(annot.base$filename, colnames(data.0))]

# adjust dys variable 
dys <- as.factor(annot.base$DYSGR_CAT)
dys <- as.character(dys)
dys[dys %in% c("CarcinomaOther", "NoData", "Unsatis")] <- "Other"
# dys[dys %in% c("Hyperplasia")] <- "Normal"
dys[dys %in% c("Hyperplasia", "Metaplasia")] <- "Normal"
# if not using all levels (2 levels: normal / abnormal)
dys[!(dys %in% c("Normal", "Other"))] <- "Dysplasia"
dys <- as.factor(dys)

mes <- data.0

d<-c()
d<-annot.base
sex <- annot.base$sex
smoke <- annot.base$smoking_status
copd <- annot.base$copd_status

var="dys"

mydata.0 <- mes

# Run limma
mod0 <- model.matrix(~ smoke + sex + copd + dys)
fit.0 <- lmFit(mydata.0, mod0)
fit.0 <- eBayes(fit.0)

res <- c()
res <- topTable(fit.0, coef="dysNormal", adjust.method="fdr", p.value=0.05, number=nrow(data.0))
# number of sig genes
cat(paste0(var, ": ", nrow(res), " / ", nrow(mydata.0), ifelse(nrow(res)>nrow(mydata.0)*0.05, " DE genes (significant; ", " DE genes (not significant; "), nrow(mydata.0)*0.05, " expected by chance)\n"))
write.csv(res, paste0("feature_selection_limma_smoke_sex_copd_dysgr_anova_seq_signal_", nrow(res), "g_lungevity.csv"))

write.table(res$ID, file=paste0("feature_selection_genes_", nrow(res), "g_lungevity.txt"), row.names=FALSE, quote=FALSE)