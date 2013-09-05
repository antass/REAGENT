setwd("~/MyDocs/BU/Kyoto/projects/aracne/data/lungevity/")
lnames <- load("~/Meta_Analysis/kyoto/aracne/data/")


all.exp.orig <- all.exp 
# Extract Ensembl genes from LUNGevity
lnames <- load("lungevity_ann_exp_82s_14308g.rda")  # lgv.ann  lgv.exp
ensembl.genes <- rownames(lgv.exp)
# write.table(ensembl.genes, file="ensembl_genes_14308g.txt", quote=FALSE, row.names=FALSE)  # use biomart to convert

# Convert Entrez IDs from BioMart into probe IDs by adding "_at"
dict.ens.ent <- read.table("biomart_ensembl_entrez_genes.txt", sep="\t", header=TRUE)
dict.ens.ent$EntrezGene.ID <- paste0(dict.ens.ent$EntrezGene.ID, "_at")

# Get list of common probes 
lgv.genes <- dict.ens.ent$EntrezGene.ID
common.probes <- Intersect(lgv.genes, rownames(all.exp), rownames(air.exp), rownames(lam.exp), rownames(grt.exp))

# Get list of all probes
union.probes <- unique(c(rownames(air.exp), rownames(lam.exp)))
union.probes <- gsub("_at", "", union.probes)
# write.table(union.probes, file="entrez_genes_19956g.txt", quote=FALSE, row.names=FALSE)  # use biomart to convert

dict.ent.ens <- read.table("biomart_entrez_ensembl_genes.txt", sep="\t", header=TRUE)
dict.ent.ens$EntrezGene.ID <- paste0(dict.ent.ens$EntrezGene.ID, "_at")

# Convert EntrezGene IDs in all 4 sets to Ensembl
# allegro
all.exp$EnsemblID <- dict.ent.ens$Ensembl.Gene.ID[ match(rownames(all.exp), dict.ent.ens$EntrezGene.ID)]
all.exp <- all.exp[, c(ncol(all.exp), 1:(ncol(all.exp)-1))]
all.exp.collapsed <- aggregate(all.exp[,-1], by=list(Gene=as.factor(all.exp$EnsemblID)), FUN=median)
save(all.exp.collapsed, file=paste0("~/Meta_Analysis/kyoto/aracne/data/all_exp_collapsed_", nrow(all.exp.collapsed), "g.rda"))

# allegro
air.exp$EnsemblID <- dict.ent.ens$Ensembl.Gene.ID[ match(rownames(air.exp), dict.ent.ens$EntrezGene.ID)]
air.exp <- air.exp[, c(ncol(air.exp), 1:(ncol(air.exp)-1))]
air.exp.collapsed <- aggregate(air.exp[,-1], by=list(Gene=as.factor(air.exp$EnsemblID)), FUN=median)
save(air.exp.collapsed, file=paste0("~/Meta_Analysis/kyoto/aracne/data/air_exp_collapsed_", nrow(air.exp.collapsed), "g.rda"))


# allegro
grt.exp$EnsemblID <- dict.ent.ens$Ensembl.Gene.ID[ match(rownames(grt.exp), dict.ent.ens$EntrezGene.ID)]
grt.exp <- grt.exp[, c(ncol(grt.exp), 1:(ncol(grt.exp)-1))]
grt.exp.collapsed <- aggregate(grt.exp[,-1], by=list(Gene=as.factor(grt.exp$EnsemblID)), FUN=median)
save(grt.exp.collapsed, file=paste0("~/Meta_Analysis/kyoto/aracne/data/grt_exp_collapsed_", nrow(grt.exp.collapsed), "g.rda"))


# allegro
all.exp$EnsemblID <- dict.ent.ens$Ensembl.Gene.ID[ match(rownames(all.exp), dict.ent.ens$EntrezGene.ID)]
all.exp <- all.exp[, c(ncol(all.exp), 1:(ncol(all.exp)-1))]
all.exp.collapsed <- aggregate(all.exp[,-1], by=list(Gene=as.factor(all.exp$EnsemblID)), FUN=median)
save(lam.exp.collapsed, file=paste0("~/Meta_Analysis/kyoto/aracne/data/lam_exp_collapsed_", nrow(lam.exp.collapsed), "g.rda"))




dict.s <- read.csv("~/MyDocs/BU/Kyoto/projects/aracne/data/feature_selction_genes_ensembl_symbol_ID_dictionary.csv")
dict.e <- read.csv("~/MyDocs/BU/Kyoto/projects/aracne/data/feature_selction_genes_ID_dictionary.csv")
dict <- merge(dict.e, dict.s, by="Ensembl.Gene.ID", all=TRUE)
dict$EntrezGene.ID <- as.factor(ifelse(!is.na(dict$EntrezGene.ID), paste0(dict$EntrezGene.ID, "_at"), NA))
dict.e$EntrezGene.ID <- as.factor(ifelse(!is.na(dict.e$EntrezGene.ID), paste0(dict.e$EntrezGene.ID, "_at"), NA))
dict.s$HGNC.symbol <- blank2na(dict.s$HGNC.symbol)



