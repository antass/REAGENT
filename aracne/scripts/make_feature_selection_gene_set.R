lnames <- load("/protected/projects/pulmarray/Airway_Metanalysis/data/processed/input/clean/Lam238_annotation_premalignancy_exprRMA_exprSCAN_clean_226.rda")  # annotation expression
lnames <- load("~/Meta_Analysis/data/processed/input/clean/Lam238_annotation_premalignancy_exprRMA_exprSCAN_clean_226.rda")  # annotation expression

# genes19741 <- rownames(expr.rma)
genes19741 <- rownames(lam.exp)

lnames <- load("/protected/projects/pulmarray/Airway_Metanalysis/data/processed/input/clean/Airway_annotation_exprRMA_exprSCAN_clean_318.rda")  # annotation expression
lnames <- load("~/Meta_Analysis/data/processed/input/clean/Airway_annotation_exprRMA_exprSCAN_clean_318.rda")  # annotation expression

# genes12146 <- rownames(expr.rma)
genes12146 <- rownames(air.exp)

commonGenes <- intersect(genes12146, genes19741)  # 11931



genes1 <- read.csv("~/Meta_Analysis/kyoto/aracne/data/feature_selection_limma_smoke_sex_cancer_anova_airway.csv")[,2]
# genes2 <- read.csv("~/Meta_Analysis/kyoto/aracne/data/feature_selection_limma_smoke_sex_cancer_anova_allegro.csv")[,2]
# genes3 <- read.csv("~/Meta_Analysis/kyoto/aracne/data/feature_selection_limma_smoke_sex_copd_dysgr_anova_greentea.csv")[,2]  # 0 genes
genes4 <- read.csv("~/Meta_Analysis/kyoto/aracne/data/feature_selection_limma_smoke_sex_copd_dysgr_anova_lam.csv")[,2]
# genes5 <- read.csv("~/kyoto/projects/aracne/data/feature_selection_limma_smoke_sex_dysgr_anova_greentea.csv")[,2]
genes6 <- read.csv("~/Meta_Analysis/kyoto/aracne/data/feature_selection_limma_smoke_sex_copd_cancer_anova_allegro.csv")[,2]



genes1 <- intersect(genes1, commonGenes)
# genes2 <- intersect(genes2, commonGenes)
# genes3 <- intersect(genes3, commonGenes)
genes4 <- intersect(genes4, commonGenes)
genes6 <- intersect(genes6, commonGenes)


kyoto.de.genes.feat.sel <- unique(c(as.character(genes1), as.character(genes4), as.character(genes6))) 
kyoto.de.genes.feat.sel.common <- Intersect(as.character(genes1), as.character(genes6), as.character(genes4)) 
kyoto.de.genes.feat.sel.atleast2 <- unique(c(intersect(genes1, genes6), intersect(genes6, genes4), intersect(genes1, genes4)))



save(kyoto.de.genes.feat.sel, file="~/Meta_Analysis/kyoto/aracne/data/Feature_selection_DE_genes_cancer_premal_union_5659.rda")
save(kyoto.de.genes.feat.sel.common, file="~/Meta_Analysis/kyoto/aracne/data/Feature_selection_DE_genes_cancer_premal_common_36.rda")
save(kyoto.de.genes.feat.sel.atleast2, file="~/Meta_Analysis/kyoto/aracne/data/Feature_selection_DE_genes_cancer_premal_atleast2_936.rda")
save(kyoto.de.genes.feat.sel, kyoto.de.genes.feat.sel.atleast2, kyoto.de.genes.feat.sel.common, 
     file="~/Meta_Analysis/kyoto/aracne/data/Feature_selection_DE_genes_all.rda")