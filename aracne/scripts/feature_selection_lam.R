library(limma, lib.loc="~/R/x86_64-unknown-linux-gnu-library/")

setwd("~/Meta_Analysis/kyoto/aracne/data/")

# ### AIRWAY

lnames <- load("/protected/projects/pulmarray/Airway_Metanalysis/data/processed/input/clean/Lam238_annotation_premalignancy_exprRMA_exprSCAN_clean_226.rda")  # annotation expression
annotation <- ann.premal
expression <- expr.rma

# annot.base <- annotation[ annotation$race %in% c("Caucasian", "African American"), ]
# annot.base$race <- droplevels(annot.base$race)
annot.base <- annotation
annot.base$copd_status <- as.character(annot.base$copd_status)
annot.base$copd_status[annot.base$copd_status=="no"] <- "no COPD"
annot.base$copd_status[annot.base$copd_status=="yes"] <- "COPD"
annot.base$copd_status <- as.factor(annot.base$copd_status)

data.0 <- expression[, as.character(annot.base$filename)]
data.0 <- data.0[, match(annot.base$filename, colnames(data.0))]

sex<-as.factor(annot.base$sex)
smoke<-as.factor(annot.base$smoking_status)
copd<-as.factor(annot.base$copd_status)

# adjust dys variable 
dys <- as.factor(annot.base$DYSGR_CAT)
dys <- as.character(dys)
dys[dys %in% c("CarcinomaOther", "NoData", "Unsatis")] <- "Other"

# if not using all levels (2 levels: normal / abnormal)
dys[!(dys %in% c("Normal", "Other"))] <- "Dysplasia"

dys <- as.factor(dys)

mes <- data.0[, !is.na(copd) & !is.na(sex) & !is.na(smoke) & dys!="Other"]
copd2 <- copd[!is.na(copd) & !is.na(sex) & !is.na(smoke) & dys!="Other"]
sex2 <- sex[!is.na(copd) & !is.na(sex) & !is.na(smoke) & dys!="Other"]
smoke2<- smoke[!is.na(copd) & !is.na(sex) & !is.na(smoke) & dys!="Other"]
dys2<- dys[!is.na(copd) & !is.na(sex) & !is.na(smoke) & dys!="Other"]


d<-c()
d<-annot.base[!is.na(copd) & !is.na(sex) & !is.na(smoke) & dys!="Other", ]

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
res <- topTable(fit.0, coef="dysNormal", adjust.method="fdr", p.value=0.05, number=nrow(data.0))
# number of sig genes
cat(paste0(var, ": ", nrow(res), " / ", nrow(mydata.0), ifelse(nrow(res)>nrow(mydata.0)*0.05, " DE genes (significant; ", " DE genes (not significant; "), nrow(mydata.0)*0.05, " expected by chance)\n"))
write.csv(res, paste0("feature_selection_limma_smoke_sex_copd_dysgr_anova_lam.csv"))
