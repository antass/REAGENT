library(limma, lib.loc="~/R/x86_64-unknown-linux-gnu-library/")

setwd("~/Meta_Analysis/kyoto/aracne/data/")

# ### AIRWAY

lnames <- load("/protected/projects/pulmarray/Airway_Metanalysis/data/processed/input/clean/GreenTea_annotation_premalignancy_exprRMA_exprSCAN_clean_61.rda")  # "ann"        "ann.premal" "expr.rma"   "expr.scan" 
annotation <- ann.premal
expression <- expr.rma
# exp <- cbind(geneID=rownames(expression), expression)
# write.table(exp, file="~/Meta_Analysis/data/processed/GreenTea_Samples_RMA_clean_61.tsv", sep="\t", row.names=FALSE)

# annot.base <- annotation[ annotation$race %in% c("Caucasian", "African American"), ]
# annot.base$race <- droplevels(annot.base$race)
annot.base <- annotation
annot.base$copd_status <- as.character(annot.base$copd_status)
annot.base$copd_status[annot.base$copd_status=="no"] <- "no COPD"
annot.base$copd_status[annot.base$copd_status=="yes"] <- "COPD"
annot.base$copd_status <- as.factor(annot.base$copd_status)
annot.base$DYSGR_CAT <- droplevels(annot.base$DYSGR_CAT)

data.0 <- expression[, as.character(annot.base$filename)]
data.0 <- data.0[, match(annot.base$filename, colnames(data.0))]

sex<-as.factor(annot.base$sex)
smoke<-as.factor(annot.base$smoking_status)
copd<-as.factor(annot.base$copd_status)

# adjust dys variable 
dys <- as.factor(annot.base$DYSGR_CAT)
dys <- as.character(dys)
dys[dys %in% c("CarcinomaOther", "NoData", "Unsatis")] <- "Other"  # will later exclude

# if not using all levels (2 levels: normal / abnormal)
dys[!(dys %in% c("Normal", "Other"))] <- "Dysplasia"

dys <- as.factor(dys)
dys <- relevel(dys, ref="Normal")

out <- !is.na(copd) & !is.na(sex) & !is.na(smoke) & dys!="Other"
mes <- data.0[, out]
copd2 <- copd[out]
sex2 <- sex[out]
smoke2<- smoke[out]
dys2<- dys[out]


d<-c()
d<-annot.base[out, ]

sex <- sex2
smoke <- droplevels(smoke2)
copd <- copd2
dys <- droplevels(dys2)

var="dys"

mydata.0 <- mes
#do Limma analysis
mod0 <- model.matrix(~ smoke + sex + copd + dys)
fit.0 <- lmFit(mydata.0, mod0, na.action=na.exclude)
fit.0 <- eBayes(fit.0)

res <- c()
#select genes by ANOVA for cancer (coef = 5 : canceryes)
res <- topTable(fit.0, coef=2, adjust.method="fdr", p.value=0.05, number=nrow(data.0))
# number of sig genes
cat(paste0(var, ": ", nrow(res), " / ", nrow(mydata.0), ifelse(nrow(res)>nrow(mydata.0)*0.05, " DE genes (significant; ", " DE genes (not significant; "), nrow(mydata.0)*0.05, " expected by chance)\n"))
write.csv(res, paste0("feature_selection_limma_smoke_sex_copd_dysgr_anova_greentea.csv"))
