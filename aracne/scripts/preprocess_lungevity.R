setwd("~/Meta_Analysis/kyoto/aracne")

# Expression
lgv.exp <- read.table("data/input/LUNGevity_log2fpkm_perc15_filter.txt", sep="\t")

lgv.baseline.samples <- as.data.frame(read.table("~/unprotected/projects/pulmseq/Premalignancy_LUNGevity/model_allbatches/base_samples.txt", sep="\t"))
lgv.samples <- gsub("Sample_", "", lgv.baseline.samples)
lgv.samples <- gsub("_accepted_hits", "", lgv.samples)
lgv.samples <- gsub("Sample_Lam_", "", lgv.samples)
lgv.samples <- gsub("Sample_Lam_", "", lgv.samples)

lgv.ann <- read.table("/unprotected/projects/pulmseq/Premalignancy_LUNGevity/annotation/LUNGevity_allsamples_curated.txt", sep="\t")

### USE LOCAL PATHS
setwd("~/MyDocs/BU/Kyoto/projects/aracne/data/lungevity/")

# Expression
lgv.exp <- read.table("LUNGevity_log2fpkm_perc15_filter.txt", sep="\t")
head(colnames(lgv.exp))

# Baseline samples
lgv.baseline.samples <- as.character(read.table("base_samples.txt", sep="\t")[,1])
# lgv.samples <- gsub("_accepted_hits", "", lgv.baseline.samples)
# lgv.samples <- gsub("\\.", "-", lgv.samples)
# lgv.samples <- gsub("Sample_", "", lgv.samples)
# lgv.samples <- gsub("Sample_Lam_", "", lgv.samples)
# lgv.samples <- gsub("Sample_Lam_", "", lgv.samples)

# Annotation
lgv.ann <- read.table("LUNGevity_allsamples_curated.txt", sep="\t", header=TRUE)
lgv.ann$filename <- paste0(lgv.ann$Seq_Name, "_accepted_hits")
lgv.ann$filename <- gsub("-", "\\.", lgv.ann$filename)

# Subset annotation and expression by baseline samples
lgv.ann.orig <- lgv.ann
lgv.ann <- lgv.ann[ lgv.ann$filename %in% lgv.baseline.samples, ]
lgv.exp <- lgv.exp[, colnames(lgv.exp) %in% lgv.baseline.samples]

## don't save to file - run lungevity_annotation_cleaning.R next


