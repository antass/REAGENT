library(limma, lib.loc="~/R/x86_64-unknown-linux-gnu-library/")

setwd("~/Meta_Analysis/kyoto/aracne/data/")
fs.genes <- read.table("features_selection_genes_4046g.txt", col.names="gene", skip=1)

# ### AIRWAY

lnames <- load("~/Meta_Analysis/data/processed/input/clean/Airway_annotation_exprRMA_exprSCAN_clean_318.rda")  # annotation expression
annotation <- ann
expression <- expr.rma[, colnames(expr.rma) %in%]

# annot.base <- annotation[ annotation$race %in% c("Caucasian", "African American"), ]
# annot.base$race <- droplevels(annot.base$race)
annot.base <- annotation
annot.base$cancer_status <- as.character(annot.base$cancer_status)
annot.base$cancer_status[annot.base$cancer_status=="no"] <- "no cancer"
annot.base$cancer_status[annot.base$cancer_status=="yes"] <- "cancer"
annot.base$cancer_status <- as.factor(annot.base$cancer_status)

data.0 <- expression[, as.character(annot.base$filename)]
data.0 <- data.0[, match(annot.base$filename, colnames(data.0))]

sex<-as.factor(annot.base$sex)
smoke<-as.factor(annot.base$smoking_status)
cancer<-as.factor(annot.base$cancer_status)

mes <- data.0[, !is.na(cancer) & !is.na(sex) & !is.na(smoke)]
cancer2 <- cancer[!is.na(cancer) & !is.na(sex) & !is.na(smoke)]
sex2 <- sex[!is.na(cancer) & !is.na(sex) & !is.na(smoke)]
smoke2<- smoke[!is.na(cancer) & !is.na(sex) & !is.na(smoke)]

d<-c()
d<-annot.base[!is.na(cancer) & !is.na(sex) & !is.na(smoke), ]

sex <- sex2
smoke <- smoke2
cancer <- cancer2

var="cancer"

mydata.0 <- mes
#do Limma analysis
mod0 <- model.matrix(~ smoke + sex + cancer)
fit.0 <- lmFit(mydata.0, mod0, na.action=na.exclude)
fit.0 <- eBayes(fit.0)

res <- c()
#select genes by ANOVA for cancer (coef = 5 : canceryes)
res <- topTable(fit.0, coef="cancerno cancer", adjust.method="fdr", p.value=0.05, number=nrow(data.0))
res.limma <- res
# number of sig genes
cat(paste0(var, ": ", nrow(res), " / ", nrow(mydata.0), ifelse(nrow(res)>nrow(mydata.0)*0.05, " DE genes (significant; ", " DE genes (not significant; "), nrow(mydata.0)*0.05, " expected by chance)\n"))
write.csv(res.limma, paste0("feature_selection_limma_smoke_sex_cancer_anova_airway.csv"))
