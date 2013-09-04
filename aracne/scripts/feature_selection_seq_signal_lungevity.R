# library(limma, lib.loc="~/R/x86_64-unknown-linux-gnu-library/")
library(limma)


setwd("~/Meta_Analysis/kyoto/aracne/data/")


lnames <- load("lungevity_ann_exp_82s_14308g.rda")  # lgv.ann  lgv.exp
dict.s <- read.csv("~/MyDocs/BU/Kyoto/projects/aracne/data/feature_selction_genes_ensembl_symbol_ID_dictionary.csv")
dict.e <- read.csv("~/MyDocs/BU/Kyoto/projects/aracne/data/feature_selction_genes_ID_dictionary.csv")
dict <- merge(dict.e, dict.s, by="Ensembl.Gene.ID", all=TRUE)
dict$EntrezGene.ID <- as.factor(ifelse(!is.na(dict$EntrezGene.ID), paste0(dict$EntrezGene.ID, "_at"), NA))
dict.e$EntrezGene.ID <- as.factor(ifelse(!is.na(dict.e$EntrezGene.ID), paste0(dict.e$EntrezGene.ID, "_at"), NA))
dict.s$HGNC.symbol <- blank2na(dict.s$HGNC.symbol)



## GET RID OF -Inf GENES (the model wouldn't run because of -Inf values)
x <- sapply(1:ncol(lgv.exp), function(i) sum(is.infinite(lgv.exp[,i])))  # can't get rid of samples, each has some -Inf
y <- sapply(1:nrow(lgv.exp), function(i) sum(is.infinite(unlist(lgv.exp[i,]))))  # lots of 0's. Great!
table(y>0)  # number of genes with -Inf
# Subset expression but excluding bad genes
lgv.exp <- lgv.exp[ which(y==0), ]

annotation <- lgv.ann
expression <- lgv.exp 

# annot.base <- annotation[ annotation$race %in% c("Caucasian", "African American"), ]
# annot.base$race <- droplevels(annot.base$race)sa
annot.base <- annotation

data.0 <- expression[, as.character(annot.base$filename)]
data.0 <- data.0[, match(annot.base$filename, colnames(data.0))]

# adjust dys variable 
dys <- as.factor(annot.base$DYSGR_CAT)
dys <- as.character(dys)
dys[dys %in% c("CarcinomaOther", "NoData", "Unsatis")] <- "Other"
dys[dys %in% c("Hyperplasia")] <- "Normal"


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
#do Limma analysis
mod0 <- model.matrix(~ smoke + sex + copd + dys)
fit.0 <- lmFit(mydata.0, mod0)
fit.0 <- eBayes(fit.0)

res <- c()
#select genes by ANOVA for cancer (coef = 5 : canceryes)
res <- topTable(fit.0, coef="dysNormal", adjust.method="fdr", p.value=0.05, number=nrow(data.0))
# number of sig genes
cat(paste0(var, ": ", nrow(res), " / ", nrow(mydata.0), ifelse(nrow(res)>nrow(mydata.0)*0.05, " DE genes (significant; ", " DE genes (not significant; "), nrow(mydata.0)*0.05, " expected by chance)\n"))
write.csv(res, paste0("feature_selection_limma_smoke_sex_copd_dysgr_anova_seq_signal_lungevity.csv"))

save(lgv.ann, lgv.exp, file="lungevity_ann_exp_82s_14126g.rda")

write.table(res$ID, file="features_selection_genes_4046g.txt", row.names=FALSE, quote=FALSE)

## Convert feature selected genes to Entrez to filter other sets
fs.genes.entrez <- dict.e$EntrezGene.ID[ match(fs.genes, dict.e$Ensembl.Gene.ID)]
fs.genes.entrez <- fs.genes.entrez[!is.na(fs.genes.entrez)]
fs.genes.entrez <- unique(fs.genes.entrez)
write.table(fs.genes.entrez, file="features_selection_genes_EntrezID_3523g.txt", quote=FALSE, row.names=FALSE)

## Take intersection of genes available to other sets and the feature selection set
fs.genes.entrez.common <- Intersect(fs.genes.entrez, rownames(air.exp), rownames(all.exp), rownames(grt.exp), rownames(lam.exp), rownames(lgv.exp))
fs.genes.ensembl.common <- dict.e$Ensembl.Gene.ID[ match(fs.genes.entrez.common, dict.e$EntrezGene.ID)]
fs.genes.ensembl.common <- unique(fs.genes.ensembl.common)

write.table(fs.genes.ensembl.common, file="feature_selection_ENSEMBL_common_2378g.txt", row.names=FALSE, quote=FALSE)
write.table(fs.genes.entrez.common, file="feature_selection_ENTREZ_common_2382g.txt", row.names=FALSE, quote=FALSE)


#