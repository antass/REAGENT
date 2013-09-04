rm(list=ls())
#load shared library path
shared.library.path <- with(R.version, sprintf("~/R_packages/R-%s.%s", major, minor));
.libPaths(c(shared.library.path, .libPaths()))

setwd("/unprotected/projects/pulmseq/Premalignancy_LUNGevity/model_allbatches/")
source("/protected/individuals/jbeane/scripts/compute_pca.R") 
library(edgeR) 
library(heatmap3) 
library(limma) 
library(sva) 
library(biomaRt)
library(GSVA)
library(org.Hs.eg.db)

sessionInfo()


#read in data and annotation
data<-read.table(file="/unprotected/projects/pulmseq/Premalignancy_LUNGevity/htseq_summary/htseq_gene_counts.txt",sep="\t",header=T,row.names=1) 
data.signal<-read.table(file="/unprotected/projects/pulmseq/Premalignancy_LUNGevity/mixture_model/htseq/htseqall_results.txt",sep="\t",header=T,row.names=1) 
annot<-read.table(file="/unprotected/projects/pulmseq/Premalignancy_LUNGevity/annotation/LUNGevity_allsamples_curated.txt",sep="\t",header=T)
lab<-paste(annot$Seq_Name,"_accepted_hits",sep="")
lab<-gsub("-",".",lab)

annot<-annot[match(colnames(data),lab),] 



gender<-read.table(file="/unprotected/projects/pulmseq/Premalignancy_LUNGevity/model/gender_genes_ensembl.txt",sep="\t") 
gender.data<-data[match(gender[,1],rownames(data)),]
gender.col<-annot$Sex 
gender.col<-gsub("F","pink",gender.col) 
gender.col<-gsub("M","blue",gender.col) 
pdf(file="heatmap_gender.pdf")
heatmap3(gender.data,ColSideColors=gender.col,col.clustering = "semisupervised",col=colorRampPalette(c("blue","white","red"), space="rgb")(255))
dev.off()

gender.gene<-gender.data[match("ENSG00000129824",rownames(gender.data)),]
males<-which(gender.gene>30)
females<-setdiff(1:ncol(gender.data),which(gender.gene>30))

gender.mismatch.ind<-c(males[which(annot$Sex[males]=="F")],females[which(annot$Sex[females]=="M")])
#rep.ind<-which(annot$Seq_Name=="Sample_8000560")
rep.ind<-which(annot$Patient=="9702197")
no.data.ind<-which(annot$DYSGR_CAT=="Unsatis")
qual.ind<-c(which(annot$Seq_Name=="Sample_101"))
#which(annot$Seq_Name=="Sample_12947"),
#which(annot$Seq_Name=="Sample_12946"))
gender.mismatch.ind
rep.ind
no.data.ind
qual.ind

bad.ind<-c(gender.mismatch.ind,rep.ind,no.data.ind,qual.ind)

bad.ind
#Subset data
annot.sub<-annot[setdiff(1:ncol(data),bad.ind),] 
data.signal.sub<-data.signal[,setdiff(1:ncol(data),bad.ind)] 
data.sub<-data[,setdiff(1:ncol(data),bad.ind)] 

#write out entire data set - without bad samples
print(dim(data.sub))
signal.dich.sub<-data.signal.sub=="signal"
signal.sum.sub<-apply(signal.dich.sub,1,sum)
perc<-0.15
sub.exp.genes<-which(signal.sum.sub>(round(perc*ncol(data.sub))))
print(length(sub.exp.genes))
data.all<-data.sub[sub.exp.genes,]
data.signal.all<-data.signal.sub[sub.exp.genes,]
data.dge.all<-DGEList(counts=data.all,genes=rownames(data.all))
data.dge.all<-calcNormFactors(data.dge.all)
cpm.data.dge.all<-cpm(data.dge.all)
write.table(cpm.data.dge.all,file="LUNGevity_log2cpm_perc15_filter.txt",sep="\t")

perc.cuts<-seq(from=0,to=1,by=0.05)
perc.res<-c()
for(p in 1:length(perc.cuts)){
	print(perc.cuts[p])
	print(length(which(signal.sum.sub>(round(perc.cuts[p]*ncol(data.sub))))))
	perc.res<-c(perc.res,length(which(signal.sum.sub>(round(perc.cuts[p]*ncol(data.sub))))))
}


#look at smoking status
smoke<-read.table(file="/unprotected/projects/pulmseq/Premalignancy_LUNGevity/model/gb_smoke_genes_ensembl.txt",sep="\t",header=T)
smoke.data<-data.sub[match(smoke[,3],rownames(data.sub)),]
smoke.col<-annot.sub$SUBGROUP
smoke.col<-gsub("EX-SMOKER","white",smoke.col)
smoke.col<-gsub("CURRENT-SMOKER","gray",smoke.col)
pdf(file="heatmap_smoke.pdf")
heatmap3(log2(smoke.data+1),ColSideColors=smoke.col,col.clustering = "semisupervised",col=colorRampPalette(c("blue","white","red"), space="rgb")(255))
dev.off()

#Chemoprevention Samples
chemo.patients<-names(which(table(annot.sub$Patient)==2))
which(annot.sub$Patient==chemo.patients)
chemo.ind<-c();
for(i in 1:length(chemo.patients)){
	chemo.ind<-c(chemo.ind,grep(chemo.patients[i],annot.sub$Patient))
}
annot.chemo<-annot.sub[chemo.ind,]
data.signal.chemo<-data.signal.sub[,chemo.ind]
data.chemo<-data.sub[,chemo.ind]
print(dim(data.chemo))
signal.dich.chemo<-data.signal.chemo=="signal"
signal.sum.chemo<-apply(signal.dich.chemo,1,sum)
perc<-0.15
chemo.exp.genes<-which(signal.sum.chemo>(round(perc*ncol(data.chemo))))
print(length(chemo.exp.genes))
data.chemo<-data.chemo[chemo.exp.genes,]
data.signal.chemo<-data.signal.chemo[chemo.exp.genes,]

#create DGEList object
data.dge.chemo<-DGEList(counts=data.chemo,genes=rownames(data.chemo))
#Scale Normalization
data.dge.chemo<-calcNormFactors(data.dge.chemo)

chemo.index<-as.factor(annot.chemo$Index)
chemo.lane<-as.factor(annot.chemo$Lane)
chemo.batch<-as.factor(annot.chemo$Batch)
chemo.sex<-as.factor(annot.chemo$Sex)
chemo.smoke<-as.factor(annot.chemo$SUBGROUP)
chemo.age<-as.numeric(annot.chemo$AGE_at_brush)
chemo.copd<-as.factor(annot.chemo$COPDdef2)
chemo.study<-as.factor(annot.chemo$StudyCategory)
chemo.dys<-as.character(annot.chemo$DYSGR_CAT)
chemo.dys<-as.factor(chemo.dys)
chemo.py<-as.numeric(annot.chemo$PACK_YEARS)
chemo.time<-annot.chemo$Base
chemo.time<-as.factor(chemo.time)
chemo.treat<-annot.chemo$Treatment
chemo.treat<-gsub("PLACEBO","NO",chemo.treat)
chemo.treat<-as.factor(chemo.treat)

#set up colors for heatmap
#set up colors for heatmap
chemo.dys.col<-chemo.dys
chemo.dys.col<-gsub("Normal","blue",chemo.dys.col)
chemo.dys.col<-gsub("Hyperplasia","cyan",chemo.dys.col)
chemo.dys.col<-gsub("Metaplasia","green",chemo.dys.col)
chemo.dys.col<-gsub("MildD","orange",chemo.dys.col)
chemo.dys.col<-gsub("ModD","pink",chemo.dys.col)
chemo.dys.col<-gsub("SevD","red",chemo.dys.col)
chemo.sex.col<-chemo.sex
chemo.sex.col<-gsub("M","blue",chemo.sex.col)
chemo.sex.col<-gsub("F","pink",chemo.sex.col)
chemo.copd.col<-chemo.copd
chemo.copd.col<-gsub(0,"white",chemo.copd.col)
chemo.copd.col<-gsub(1,"black",chemo.copd.col)
chemo.study.col<-chemo.study
chemo.study.col<-gsub("GTE","green",chemo.study.col)
chemo.study.col<-gsub("MYOIN","purple",chemo.study.col)
chemo.study.col<-gsub("BDSUL","yellow",chemo.study.col)
chemo.study.col<-gsub("NLSS","cyan",chemo.study.col)
chemo.smoke.col<-chemo.smoke
chemo.smoke.col<-gsub("CURRENT-SMOKER","red",chemo.smoke.col)
chemo.smoke.col<-gsub("EX-SMOKER","green",chemo.smoke.col)
chemo.time.col<-chemo.time
chemo.time.col<-gsub("post","gray",chemo.time.col)
chemo.time.col<-gsub("base","black",chemo.time.col)
chemo.treat.col<-chemo.treat
chemo.treat.col<-gsub("GTE","green",chemo.treat.col)
chemo.treat.col<-gsub("MYO","purple",chemo.treat.col)
chemo.treat.col<-gsub("NO","white",chemo.treat.col)

#Baseline analysis -- find changes associated wtih grade
annot.base<-annot.sub[which(annot.sub$Base=="base"),]
data.signal.base<-data.signal.sub[,which(annot.sub$Base=="base")]
data.base<-data.sub[,which(annot.sub$Base=="base")]
print(dim(data.base))
signal.dich.base<-data.signal.base=="signal"
signal.sum.base<-apply(signal.dich.base,1,sum)
perc<-0.15
base.exp.genes<-which(signal.sum.base>(round(perc*ncol(data.base))))
print(length(base.exp.genes))
data.base<-data.base[base.exp.genes,]
data.signal.base<-data.signal.base[base.exp.genes,]

#create DGEList object
data.dge<-DGEList(counts=data.base,genes=rownames(data.base))
#Scale Normalization
data.dge<-calcNormFactors(data.dge)

index<-as.factor(annot.base$Index)
lane<-as.factor(annot.base$Lane)
batch<-as.factor(annot.base$Batch)
sex<-as.factor(annot.base$Sex)
smoke<-as.factor(annot.base$SUBGROUP)
age<-as.numeric(annot.base$AGE_at_brush)
copd<-as.factor(annot.base$COPDdef2)
study<-as.factor(annot.base$StudyCategory)
dys<-as.character(annot.base$DYSGR_CAT)
dys<-as.factor(dys)
py<-as.numeric(annot.base$PACK_YEARS)


#set up colors for heatmap
dys.col<-dys
dys.col<-gsub("Normal","blue",dys.col)
dys.col<-gsub("Hyperplasia","cyan",dys.col)
dys.col<-gsub("Metaplasia","green",dys.col)
dys.col<-gsub("MildD","orange",dys.col)
dys.col<-gsub("ModD","pink",dys.col)
dys.col<-gsub("SevD","red",dys.col)
sex.col<-sex
sex.col<-gsub("M","blue",sex.col)
sex.col<-gsub("F","pink",sex.col)
copd.col<-copd
copd.col<-gsub(0,"white",copd.col)
copd.col<-gsub(1,"black",copd.col)
study.col<-study
study.col<-gsub("GTE","green",study.col)
study.col<-gsub("MYOIN","purple",study.col)
study.col<-gsub("BDSUL","yellow",study.col)
study.col<-gsub("NLSS","cyan",study.col)
smoke.col<-smoke
smoke.col<-gsub("CURRENT-SMOKER","red",smoke.col)
smoke.col<-gsub("EX-SMOKER","green",smoke.col)

mod0<-model.matrix(~sex+smoke+copd+dys)
data.0<-voom(data.dge,mod0,plot=TRUE)
dev.off()
data.v<-data.0
pdf(file="basesamples_MDS.pdf")
plotMDS(data.v,gene.selection="common",labels=rep(1,length(ncol(data.base))))
dev.off()
pdf(file="basesamples_MDS_batch.pdf")
plotMDS(data.v,gene.selection="common",labels=batch)
dev.off()
pdf(file="basesamples_MDS_index.pdf")
plotMDS(data.v,gene.selection="common",labels=index)
dev.off()
pdf(file="basesamples_MDS_lane.pdf")
plotMDS(data.v,gene.selection="common",labels=lane)
dev.off()
pdf(file="basesamples_MDS_smoke.pdf")
plotMDS(data.v,gene.selection="common",labels=smoke)
dev.off()
pdf(file="basesamples_MDS_study.pdf")
plotMDS(data.v,gene.selection="common",labels=study)
dev.off()
pdf(file="basesamples_MDS_copd.pdf")
plotMDS(data.v,gene.selection="common",labels=copd)
dev.off()
pdf(file="basesamples_MDS_sex.pdf")
plotMDS(data.v,gene.selection="common",labels=sex)
dev.off()
pdf(file="basesamples_MDS_dys.pdf")
plotMDS(data.v,gene.selection="common",labels=dys)
dev.off()
pdf(file="basesamples_MDS_size.pdf")
plotMDS(data.v,gene.selection="common",labels=trunc(data.dge$samples[,2]/1000000))
dev.off()
pdf(file="basesamples_MDS_rin.pdf")
plotMDS(data.v,gene.selection="common",labels=annot.base$Rin)
dev.off()
pdf(file="basesamples_top50.pdf")
plotMDS(data.v,top=50,labels=batch,gene.selection="common")
dev.off()
pdf(file="basesamples_MDS_py.pdf")
plotMDS(data.v,gene.selection="common",labels=annot.base$PACK_YEARS)
dev.off()

pdf(file="basesamples_MDS_top5000.pdf")
plotMDS(data.v,top=5000,gene.selection="common",labels=rep(1,length(ncol(data.base))))
dev.off()
pdf(file="basesamples_MDS_batch_top5000.pdf")
plotMDS(data.v,top=5000,gene.selection="common",labels=batch)
dev.off()
pdf(file="basesamples_MDS_index_top5000.pdf")
plotMDS(data.v,top=5000,gene.selection="common",labels=index)
dev.off()
pdf(file="basesamples_MDS_lane_top5000.pdf")
plotMDS(data.v,top=5000,gene.selection="common",labels=lane)
dev.off()
pdf(file="basesamples_MDS_smoke_top5000.pdf")
plotMDS(data.v,top=5000,gene.selection="common",labels=smoke)
dev.off()
pdf(file="basesamples_MDS_study_top5000.pdf")
plotMDS(data.v,top=5000,gene.selection="common",labels=study)
dev.off()
pdf(file="basesamples_MDS_copd_top5000.pdf")
plotMDS(data.v,top=5000,gene.selection="common",labels=copd)
dev.off()
pdf(file="basesamples_MDS_sex_top5000.pdf")
plotMDS(data.v,top=5000,gene.selection="common",labels=sex)
dev.off()
pdf(file="basesamples_MDS_dys_top5000.pdf")
plotMDS(data.v,top=5000,gene.selection="common",labels=dys)
dev.off()
pdf(file="basesamples_MDS_size_top5000.pdf")
plotMDS(data.v,top=5000,gene.selection="common",labels=trunc(data.dge$samples[,2]/1000000))
dev.off()
pdf(file="basesamples_MDS_rin_top5000.pdf")
plotMDS(data.v,top=5000,gene.selection="common",labels=annot.base$Rin)
dev.off()
pdf(file="basesamples_MDS_py_top5000.pdf")
plotMDS(data.v,top=5000,gene.selection="common",labels=annot.base$PACK_YEARS)
dev.off()

pdf(file="basesamples_MDS_allgenes_sex.pdf")
plotMDS(data.v,gene.selection="common",labels=sex,top=13757)
dev.off()
#pdf(file="baseamples_MDS_nontrans_logcountdata.pdf")
#mds.logcount<-plotMDS(log2(data.base+1),gene.selection="common",labels=sex,top=13757)
#dev.off()
pdf(file="basesamples_PCA_nontrans_datadge.pdf")
z.1<-c()
z.1<-(data.dge$counts-apply(data.dge$counts,1,mean))/apply(data.dge$counts,1,sd)
pca.1<-c()
pca.1<-compute_pca(z.1)
plot(pca.1[,1],pca.1[,2],pch=20,col="blue")
points(pca.1[which(sex=="F"),1],pca.1[which(sex=="F"),2],pch=20,col="pink")
dev.off()
pdf(file="baseamples_PCA_nontrans_logcountdata.pdf")
z.1<-c()
z.1<-(log2(data.base+1)-apply(log2(data.base+1),1,mean))/apply(log2(data.base+1),1,sd)
pca.1<-c()
pca.1<-compute_pca(z.1)
pca.logcount<-prcomp(as.matrix(z.1),center=F)$rotation
plot(pca.1[,1],pca.1[,2],pch=20,col="blue")
points(pca.1[which(sex=="F"),1],pca.1[which(sex=="F"),2],pch=20,col="pink")
dev.off()

dys.cont<-dys
dys.cont<-gsub("Normal",0,dys.cont)
dys.cont<-gsub("Hyperplasia",1,dys.cont)
dys.cont<-gsub("Metaplasia",2,dys.cont)
dys.cont<-gsub("MildD",3,dys.cont)
dys.cont<-gsub("ModD",4,dys.cont)
dys.cont<-gsub("SevD",5,dys.cont)
dys.cont<-as.numeric(dys.cont)
dys<-as.factor(dys.cont)

pval<-c(0.00001,0.0001,0.001,0.01)

#do Limma analysis
mod0<-model.matrix(~sex+smoke+copd+dys)
data.0<-voom(data.dge,mod0,plot=FALSE)
fit.0<-lmFit(data.0,mod0)
fit.0<-eBayes(fit.0)
res<-c()
#select genes by ANOVA for dys
res<-topTable(fit.0,coef=5:9,adjust.method="fdr",p.value=0.05,number=nrow(data.0))
print(dim(res))
write.table(res, "base_dysfact_wsex_smoke_copd_fdr05.txt", sep="\t")
res1<-res

model.sig.genes<-c()
for(i in 1:length(pval)){model.sig.genes<-c(model.sig.genes,length(which(res$P.Value<pval[i])))}


#do Limma analysis
mod0<-model.matrix(~sex+smoke+copd+dys.cont)
data.0<-voom(data.dge,mod0,plot=FALSE)
fit.0<-lmFit(data.0,mod0)
fit.0<-eBayes(fit.0)
res<-c()
#select genes by ANOVA for dys
res<-topTable(fit.0,coef=5,adjust.method="fdr",p.value=0.05,number=nrow(data.0))
print(dim(res))
write.table(res, "base_dyscont_wsex_smoke_copd_fdr05.txt", sep="\t")
res2<-res
for(i in 1:length(pval)){model.sig.genes<-c(model.sig.genes,length(which(res$P.Value<pval[i])))}


length(intersect(res1$genes,res2$genes))
#647


#do Limma analysis
mod0<-model.matrix(~sex+smoke+copd+study+dys)
data.0<-voom(data.dge,mod0,plot=FALSE)
fit.0<-lmFit(data.0,mod0)
fit.0<-eBayes(fit.0)
res<-c()
#select genes by ANOVA for dys
res<-topTable(fit.0,coef=8:12,adjust.method="fdr",p.value=0.05,number=nrow(data.0))
print(dim(res))
write.table(res, "base_dysfact_wsex_smoke_copd_study_fdr05.txt", sep="\t")
for(i in 1:length(pval)){model.sig.genes<-c(model.sig.genes,length(which(res$P.Value<pval[i])))}


#do Limma analysis
mod0<-model.matrix(~sex+smoke+study+dys)
data.0<-voom(data.dge,mod0,plot=FALSE)
fit.0<-lmFit(data.0,mod0)
fit.0<-eBayes(fit.0)
res<-c()
#select genes by ANOVA for dys
res<-topTable(fit.0,coef=7:11,adjust.method="fdr",p.value=0.05,number=nrow(data.0))
print(dim(res))
write.table(res, "base_dysfact_wsex_smoke_study_fdr05.txt", sep="\t")
for(i in 1:length(pval)){model.sig.genes<-c(model.sig.genes,length(which(res$P.Value<pval[i])))}



#do Limma analysis
mod0<-model.matrix(~sex+smoke+dys)
data.0<-voom(data.dge,mod0,plot=FALSE)
fit.0<-lmFit(data.0,mod0)
fit.0<-eBayes(fit.0)
res<-c()
#select genes by ANOVA for dys
res<-topTable(fit.0,coef=4:8,adjust.method="fdr",p.value=0.05,number=nrow(data.0))
print(dim(res))
write.table(res, "base_dysfact_wsex_smoke_fdr05.txt", sep="\t")
for(i in 1:length(pval)){model.sig.genes<-c(model.sig.genes,length(which(res$P.Value<pval[i])))}


#do Limma analysis
mod0<-model.matrix(~sex+smoke+py+dys)
data.0<-voom(data.dge,mod0,plot=FALSE)
fit.0<-lmFit(data.0,mod0)
fit.0<-eBayes(fit.0)
res<-c()
#select genes by ANOVA for dys
res<-topTable(fit.0,coef=5:9,adjust.method="fdr",p.value=0.05,number=nrow(data.0))
print(dim(res))
write.table(res, "base_dysfact_wsex_smoke_py_fdr05.txt", sep="\t")
for(i in 1:length(pval)){model.sig.genes<-c(model.sig.genes,length(which(res$P.Value<pval[i])))}



#do Limma analysis
mod0<-model.matrix(~sex+smoke+age+dys)
data.0<-voom(data.dge,mod0,plot=FALSE)
fit.0<-lmFit(data.0,mod0)
fit.0<-eBayes(fit.0)
res<-c()
#select genes by ANOVA for dys
res<-topTable(fit.0,coef=5:9,adjust.method="fdr",p.value=0.05,number=nrow(data.0))
print(dim(res))
write.table(res, "base_dysfact_wsex_smoke_age_fdr05.txt", sep="\t")
for(i in 1:length(pval)){model.sig.genes<-c(model.sig.genes,length(which(res$P.Value<pval[i])))}



#do Limma analysis
mod0<-model.matrix(~sex+smoke+batch+dys)
data.0<-voom(data.dge,mod0,plot=FALSE)
fit.0<-lmFit(data.0,mod0)
fit.0<-eBayes(fit.0)
res<-c()
#select genes by ANOVA for dys
res<-topTable(fit.0,coef=7:11,adjust.method="fdr",p.value=0.05,number=nrow(data.0))
print(dim(res))
write.table(res, "base_dysfact_wsex_smoke_batch_fdr05.txt", sep="\t")
for(i in 1:length(pval)){model.sig.genes<-c(model.sig.genes,length(which(res$P.Value<pval[i])))}


#do Limma analysis
mod0<-model.matrix(~sex+smoke+lane+dys)
data.0<-voom(data.dge,mod0,plot=FALSE)
fit.0<-lmFit(data.0,mod0)
fit.0<-eBayes(fit.0)
res<-c()
#select genes by ANOVA for dys
res<-topTable(fit.0,coef=11:15,adjust.method="fdr",p.value=0.05,number=nrow(data.0))
print(dim(res))
write.table(res, "base_dysfact_wsex_smoke_lane_fdr05.txt", sep="\t")
for(i in 1:length(pval)){model.sig.genes<-c(model.sig.genes,length(which(res$P.Value<pval[i])))}


#do Limma analysis
mod0<-model.matrix(~sex+smoke+index+dys)
data.0<-voom(data.dge,mod0,plot=FALSE)
fit.0<-lmFit(data.0,mod0)
fit.0<-eBayes(fit.0)
res<-c()
#select genes by ANOVA for dys
res<-topTable(fit.0,coef=7:11,adjust.method="fdr",p.value=0.05,number=nrow(data.0))
print(dim(res))
write.table(res, "base_dysfact_wsex_smoke_index_fdr05.txt", sep="\t")
for(i in 1:length(pval)){model.sig.genes<-c(model.sig.genes,length(which(res$P.Value<pval[i])))}

write.table(model.sig.genes,file="models_sig_genes.txt",sep="\t")

############
dys.cont2<-as.character(annot.base$DYSGR_CAT)
dys.cont2<-gsub("Normal",0,dys.cont2)
dys.cont2<-gsub("Hyperplasia",0,dys.cont2)
dys.cont2<-gsub("Metaplasia",0,dys.cont2)
dys.cont2<-gsub("MildD",1,dys.cont2)
dys.cont2<-gsub("ModD",2,dys.cont2)
dys.cont2<-gsub("SevD",3,dys.cont2)
dys.cont2<-as.numeric(dys.cont2)
dys2<-as.factor(dys.cont2)
mod0<-model.matrix(~sex+smoke+copd+dys2)
data.0<-voom(data.dge,mod0,plot=FALSE)
fit.0<-lmFit(data.0,mod0)
fit.0<-eBayes(fit.0)
res<-c()
res<-topTable(fit.0,coef=5:7,adjust.method="fdr",p.value=0.05,number=nrow(data.0))
write.table(res,file="dysfact_wsex_smoke_copd_l4_fdr05.txt",sep="\t")
res<-topTable(fit.0,coef=5:7,adjust.method="fdr",p.value=1,number=nrow(data.0))
write.table(res,file="dysfact_wsex_smoke_copd_l4.txt",sep="\t")
res.l4<-res

dys.cont2<-as.character(annot.base$DYSGR_CAT)
dys.cont2<-gsub("Normal",0,dys.cont2)
dys.cont2<-gsub("Hyperplasia",0,dys.cont2)
dys.cont2<-gsub("Metaplasia",0,dys.cont2)
dys.cont2<-gsub("MildD",1,dys.cont2)
dys.cont2<-gsub("ModD",1,dys.cont2)
dys.cont2<-gsub("SevD",1,dys.cont2)
dys.cont2<-as.numeric(dys.cont2)
dys2<-as.factor(dys.cont2)
mod0<-model.matrix(~sex+smoke+copd+dys2)
data.0<-voom(data.dge,mod0,plot=FALSE)
fit.0<-lmFit(data.0,mod0)
fit.0<-eBayes(fit.0)
res<-c()
res<-topTable(fit.0,coef=5,adjust.method="fdr",p.value=0.05,number=nrow(data.0))
write.table(res,file="dysfact_wsex_smoke_copd_l2_fdr05.txt",sep="\t")
res<-topTable(fit.0,coef=5,adjust.method="fdr",p.value=1,number=nrow(data.0))
write.table(res,file="dysfact_wsex_smoke_copd_l2.txt",sep="\t")
res.l2<-res
res<-c();
res<-topTable(fit.0,coef=5,adjust.method="fdr",p.value=1,number=nrow(data.0),lfc=log2(1.5))
write.table(res,file="dysfact_wsex_smoke_copd_l2_fc1.5.txt",sep="\t")
res.l2.fc<-res
res.l2.sign<-res.l2$t
res.l2.sign[which(res.l2.sign<0)]<-(-1)
res.l2.sign[which(res.l2.sign>=0)]<-(1)
table(res.l2.sign[intersect(which(abs(res.l2$logFC)>log2(1.5)),which(res.l2$adj.P.Val<0.01))]) 
table(res.l2.sign[intersect(which(abs(res.l2$logFC)>log2(1)),which(res.l2$adj.P.Val<0.01))])
table(res.l2.sign[intersect(which(abs(res.l2$logFC)>log2(1.5)),which(res.l2$adj.P.Val<0.05))])
table(res.l2.sign[intersect(which(abs(res.l2$logFC)>log2(1)),which(res.l2$adj.P.Val<0.05))])
table(res.l2.sign[intersect(which(abs(res.l2$logFC)>log2(1.5)),which(res.l2$adj.P.Val<0.01))])
res.test<-topTable(fit.0,coef=5,adjust.method="BH",p.value=1,number=nrow(data.0))
res.test2<-topTable(fit.0,coef=5,adjust.method="fdr",p.value=1,number=nrow(data.0))
dys.cont2<-as.character(annot.base$DYSGR_CAT)
dys.cont2<-gsub("Normal",0,dys.cont2)
dys.cont2<-gsub("Hyperplasia",0,dys.cont2)
dys.cont2<-gsub("Metaplasia",1,dys.cont2)
dys.cont2<-gsub("MildD",2,dys.cont2)
dys.cont2<-gsub("ModD",2,dys.cont2)
dys.cont2<-gsub("SevD",2,dys.cont2)
dys.cont2<-as.numeric(dys.cont2)
dys2<-as.factor(dys.cont2)
mod0<-model.matrix(~sex+smoke+copd+dys2)
data.0<-voom(data.dge,mod0,plot=FALSE)
fit.0<-lmFit(data.0,mod0)
fit.0<-eBayes(fit.0)
res<-c()
res<-topTable(fit.0,coef=5:6,adjust.method="fdr",p.value=0.05,number=nrow(data.0))
write.table(res,file="dysfact_wsex_smoke_copd_l3wmeta_fdr05.txt",sep="\t")
res<-topTable(fit.0,coef=5:6,adjust.method="fdr",p.value=1,number=nrow(data.0))
write.table(res,file="dysfact_wsex_smoke_copd_l3wmeta.txt",sep="\t")
res.l3<-res

#############
#Final model for now
mod0<-model.matrix(~sex+smoke+copd+dys)
data.0<-voom(data.dge,mod0,plot=FALSE)
fit.0<-lmFit(data.0,mod0)
fit.0<-eBayes(fit.0)
res<-c()

#select genes by ANOVA for dys
res.fc<-topTable(fit.0,coef=5:9,adjust.method="fdr",p.value=0.05,number=nrow(data.0),lfc=log2(1.5))
write.table(res.fc,file="dysfact_wsex_smoke_copd_fdr05_fc1pt5.txt",sep="\t")

res<-topTable(fit.0,coef=5:9,adjust.method="fdr",p.value=0.05,number=nrow(data.0))
write.table(res,file="dysfact_wsex_smoke_copd_fdr05.txt",sep="\t")

res<-topTable(fit.0,coef=5:9,adjust.method="fdr",p.value=1,number=nrow(data.0))
write.table(res,file="dysfact_wsex_smoke_copd.txt",sep="\t")

model.sig.genes<-c()
for(i in 1:length(pval)){model.sig.genes<-c(model.sig.genes,length(which(res$P.Value<pval[i])))}

#make results in heatmaps one with foldchange cutoff
#res<-res.fc

#compute residuals for data
mod.r<-model.matrix(~sex+smoke+copd)
data.r<-voom(data.dge,mod.r,plot=FALSE)
fit.r<-lmFit(data.r,mod.r)
fit.r<-eBayes(fit.r)
resid<-residuals(fit.r,data.r)

pdf(file="dys_wsmoke_sex_copd_volcano.pdf")
volcanoplot(fit.0,coef=5:9,highlight=149)
dev.off()

cpm.data.dge<-cpm(data.dge)
cpm.chem.dge<-cpm(data.dge.chemo)

chemo.column.matrix<-matrix(rbind(chemo.dys.col,chemo.time.col,chemo.treat.col,chemo.smoke.col,chemo.sex.col,chemo.copd.col,chemo.study.col),nrow=7,ncol=length(chemo.smoke.col))
column.matrix<-matrix(rbind(dys.col,smoke.col,sex.col,copd.col,study.col),nrow=5,ncol=length(smoke.col))
pdf(file="dys_wsmoke_sex_copd_fdr01_log2cpm.pdf")
sig.genes<-res$genes[which(res$adj.P.Val<0.01)]
data.heatmap<-cpm.data.dge[match(sig.genes,rownames(data.dge)),]
heatmap3(data.heatmap,ColSideColors=as.matrix(t(column.matrix)),col.clustering = "unsupervised",col=colorRampPalette(c("blue","white","red"), space="rgb")(255))
dev.off()
pdf(file="dys_wsmoke_sex_copd_fdr01_voomE.pdf")
sig.genes<-res$genes[which(res$adj.P.Val<0.01)]
data.heatmap<-data.0$E[match(sig.genes,rownames(data.0)),]
heatmap3(data.heatmap,ColSideColors=as.matrix(t(column.matrix)),col.clustering = "unsupervised",col=colorRampPalette(c("blue","white","red"), space="rgb")(255))
dev.off()
pdf(file="dys_wsmoke_sex_copd_fdr01_resid.pdf")
sig.genes<-res$genes[which(res$adj.P.Val<0.01)]
data.heatmap<-resid[match(sig.genes,rownames(resid)),]
heatmap3(data.heatmap,ColSideColors=as.matrix(t(column.matrix)),col.clustering = "unsupervised",col=colorRampPalette(c("blue","white","red"), space="rgb")(255))
dev.off()
pdf(file="chemosamples_dys_wsmoke_sex_copd_fdr01_log2cpm.pdf")
sig.genes<-res$genes[which(res$adj.P.Val<0.01)]
data.heatmap<-cpm.chem.dge[match(sig.genes,rownames(data.dge.chemo))[which(is.na(match(sig.genes,rownames(data.dge.chemo)))==FALSE)],]
heatmap3(data.heatmap,ColSideColors=as.matrix(t(chemo.column.matrix)),col.clustering = "unsupervised",col=colorRampPalette(c("blue","white","red"), space="rgb")(255))
dev.off()

pdf(file="dys_wsmoke_sex_copd_semi_fdr01_log2cpm.pdf")
sig.genes<-res$genes[which(res$adj.P.Val<0.01)]
length(sig.genes)
data.heatmap<-cpm.data.dge[match(sig.genes,rownames(data.dge)),]
heatmap3(data.heatmap,ColSideColors=as.matrix(t(column.matrix)),col.clustering = "semisupervised",col=colorRampPalette(c("blue","white","red"), space="rgb")(255))
dev.off()

pdf(file="dys_wsmoke_sex_copd_semi_fdr01_voomE.pdf")
sig.genes<-res$genes[which(res$adj.P.Val<0.01)]
data.heatmap<-data.0$E[match(sig.genes,rownames(data.0)),]
heatmap3(data.heatmap,ColSideColors=as.matrix(t(column.matrix)),col.clustering = "semisupervised",col=colorRampPalette(c("blue","white","red"), space="rgb")(255))
dev.off()
pdf(file="dys_wsmoke_sex_copd_semi_fdr01_resid.pdf")
sig.genes<-res$genes[which(res$adj.P.Val<0.01)]
data.heatmap<-resid[match(sig.genes,rownames(resid)),]
heatmap3(data.heatmap,ColSideColors=as.matrix(t(column.matrix)),col.clustering = "semisupervised",col=colorRampPalette(c("blue","white","red"), space="rgb")(255))
dev.off()
pdf(file="chemosamples_dys_wsmoke_sex_copd_semi_fdr01_log2cpm.pdf")
sig.genes<-res$genes[which(res$adj.P.Val<0.01)]
data.heatmap<-cpm.chem.dge[match(sig.genes,rownames(data.dge.chemo))[which(is.na(match(sig.genes,rownames(data.dge.chemo)))==FALSE)],]
heatmap3(data.heatmap,ColSideColors=as.matrix(t(chemo.column.matrix)),col.clustering = "semisupervised",col=colorRampPalette(c("blue","white","red"), space="rgb")(255))
dev.off()

pdf(file="dys_wsmoke_sex_copd_semi_fdr05_log2cpm.pdf")
sig.genes<-res$genes[which(res$adj.P.Val<0.05)]
length(sig.genes)
data.heatmap<-cpm.data.dge[match(sig.genes,rownames(data.dge)),]
heatmap3(data.heatmap,ColSideColors=as.matrix(t(column.matrix)),col.clustering = "semisupervised",col=colorRampPalette(c("blue","white","red"), space="rgb")(255))
dev.off()
pdf(file="dys_wsmoke_sex_copd_semi_fdr05_voomE.pdf")
sig.genes<-res$genes[which(res$adj.P.Val<0.05)]
data.heatmap<-data.0$E[match(sig.genes,rownames(data.0)),]
heatmap3(data.heatmap,ColSideColors=as.matrix(t(column.matrix)),col.clustering = "semisupervised",col=colorRampPalette(c("blue","white","red"), space="rgb")(255))
dev.off()
pdf(file="dys_wsmoke_sex_copd_semi_fdr05_resid.pdf")
sig.genes<-res$genes[which(res$adj.P.Val<0.05)]
data.heatmap<-resid[match(sig.genes,rownames(resid)),]
heatmap3(data.heatmap,ColSideColors=as.matrix(t(column.matrix)),col.clustering = "semisupervised",col=colorRampPalette(c("blue","white","red"), space="rgb")(255))
dev.off()
pdf(file="chemosamples_dys_wsmoke_sex_copd_semi_fdr05_log2cpm.pdf")
sig.genes<-res$genes[which(res$adj.P.Val<0.05)]
data.heatmap<-cpm.chem.dge[match(sig.genes,rownames(data.dge.chemo))[which(is.na(match(sig.genes,rownames(data.dge.chemo)))==FALSE)],]
heatmap3(data.heatmap,ColSideColors=as.matrix(t(chemo.column.matrix)),col.clustering = "semisupervised",col=colorRampPalette(c("blue","white","red"), space="rgb")(255))
dev.off()

pdf(file="dys_wsmoke_sex_copd_fdr05_log2cpm.pdf")
sig.genes<-res$genes[which(res$adj.P.Val<0.05)]
length(sig.genes)
data.heatmap<-cpm.data.dge[match(sig.genes,rownames(data.dge)),]
heatmap3(data.heatmap,ColSideColors=as.matrix(t(column.matrix)),col.clustering = "unsupervised",col=colorRampPalette(c("blue","white","red"), space="rgb")(255))
dev.off()
pdf(file="dys_wsmoke_sex_copd_fdr05_voomE.pdf")
sig.genes<-res$genes[which(res$adj.P.Val<0.05)]
data.heatmap<-data.0$E[match(sig.genes,rownames(data.0)),]
heatmap3(data.heatmap,ColSideColors=as.matrix(t(column.matrix)),col.clustering = "unsupervised",col=colorRampPalette(c("blue","white","red"), space="rgb")(255))
dev.off()
pdf(file="dys_wsmoke_sex_copd_fdr05_resid.pdf")
sig.genes<-res$genes[which(res$adj.P.Val<0.05)]
data.heatmap<-resid[match(sig.genes,rownames(resid)),]
heatmap3(data.heatmap,ColSideColors=as.matrix(t(column.matrix)),col.clustering = "unsupervised",col=colorRampPalette(c("blue","white","red"), space="rgb")(255))
dev.off()
pdf(file="chemosamples_dys_wsmoke_sex_copd_fdr05_log2cpm.pdf")
sig.genes<-res$genes[which(res$adj.P.Val<0.05)]
data.heatmap<-cpm.chem.dge[match(sig.genes,rownames(data.dge.chemo))[which(is.na(match(sig.genes,rownames(data.dge.chemo)))==FALSE)],]
heatmap3(data.heatmap,ColSideColors=as.matrix(t(chemo.column.matrix)),col.clustering = "unsupervised",col=colorRampPalette(c("blue","white","red"), space="rgb")(255))
dev.off()

pdf(file="dys_wsmoke_sex_copd_semi_l2_fdr05_log2cpm.pdf")
sig.genes<-res.l2$genes[which(res.l2$adj.P.Val<0.05)]
length(sig.genes)
data.heatmap<-cpm.data.dge[match(sig.genes,rownames(data.dge)),]
heatmap3(data.heatmap,ColSideColors=as.matrix(t(column.matrix)),col.clustering = "semisupervised",col=colorRampPalette(c("blue","white","red"), space="rgb")(255))
dev.off()
pdf(file="dys_wsmoke_sex_copd_semi_l2_fdr01_log2cpm.pdf")
sig.genes<-res.l2$genes[which(res.l2$adj.P.Val<0.01)]
length(sig.genes)
data.heatmap<-cpm.data.dge[match(sig.genes,rownames(data.dge)),]
heatmap3(data.heatmap,ColSideColors=as.matrix(t(column.matrix)),col.clustering = "semisupervised",col=colorRampPalette(c("blue","white","red"), space="rgb")(255))
dev.off()
pdf(file="dys_wsmoke_sex_copd_semi_l2_fdr001_log2cpm.pdf")
sig.genes<-res.l2$genes[which(res.l2$adj.P.Val<0.001)]
length(sig.genes)
data.heatmap<-cpm.data.dge[match(sig.genes,rownames(data.dge)),]
heatmap3(data.heatmap,ColSideColors=as.matrix(t(column.matrix)),col.clustering = "semisupervised",col=colorRampPalette(c("blue","white","red"), space="rgb")(255))
dev.off()

pdf(file="dys_wsmoke_sex_copd_semi_l3_fdr05_log2cpm.pdf")
sig.genes<-res.l3$genes[which(res.l3$adj.P.Val<0.05)]
length(sig.genes)
data.heatmap<-cpm.data.dge[match(sig.genes,rownames(data.dge)),]
heatmap3(data.heatmap,ColSideColors=as.matrix(t(column.matrix)),col.clustering = "semisupervised",col=colorRampPalette(c("blue","white","red"), space="rgb")(255))
dev.off()
pdf(file="dys_wsmoke_sex_copd_semi_l4_fdr05_log2cpm.pdf")
sig.genes<-res.l4$genes[which(res.l4$adj.P.Val<0.05)]
length(sig.genes)
data.heatmap<-cpm.data.dge[match(sig.genes,rownames(data.dge)),]
heatmap3(data.heatmap,ColSideColors=as.matrix(t(column.matrix)),col.clustering = "semisupervised",col=colorRampPalette(c("blue","white","red"), space="rgb")(255))
dev.off()

####GSVA

#read in ensembl to entrez gene mapping
entrez<-read.table(file="/unprotected/projects/pulmseq/Premalignancy_LUNGevity/model/ensembl_toentrez.txt",sep="\t",header=T)
data.comp.exp<-c();
data.comp.exp<-data.dge$counts
entrez.id<-entrez[match(rownames(data.comp.exp),entrez[,1]),2]
d.entrez<-data.comp.exp[which(is.na(entrez.id)==FALSE),]
entrez.id<-entrez.id[which(is.na(entrez.id)==FALSE)]
red.id<-names(which(table(entrez.id[which(is.na(entrez.id)==FALSE)])>1))
nonred.id<-setdiff(entrez.id,red.id)
ind.final<-c()
for(i in 1:length(red.id)){
        ind<-which(entrez.id==red.id[i])
        mean.ind<-apply(d.entrez[ind,],1,mean)
        good.ind<-which(mean.ind==max(mean.ind))
        ind.final<-c(ind.final,ind[good.ind])
}
all.ind<-c(match(nonred.id,entrez.id),ind.final)
lab<-c(nonred.id,entrez.id[ind.final])
d.final<-d.entrez[all.ind,]
rownames(d.final)<-lab

pheno<-cbind(sex,copd,smoke)
rownames(pheno)<-colnames(d.final)
pheno.var<-c("sex","copd","smoke")
pheno.desc<-c("sex","copd","smoke")
pheno.df<-data.frame(pheno.var,pheno.desc)
pd <- new("AnnotatedDataFrame", data=data.frame(pheno),varMetadata=pheno.df)

e.set<-new("ExpressionSet",exprs=as.matrix(d.final),phenoData=pd,annotation="org.Hs.eg.db")
#####read in GMT files########
gene.set<-getGmt("/unprotected/projects/pulmseq/Premalignancy_LUNGevity/model/c2.cp.v3.1.entrez.gmt",geneIdType=EntrezIdentifier())
#GSVA values
gsva.val <- gsva(e.set,gene.set,min.sz=10, max.sz=500,verbose=TRUE,rnaseq=TRUE)$es.obs
write.table(gsva.val,file="htseq_GSVAvals.txt",sep="\t")

#Final model for now
mod0<-model.matrix(~sex+smoke+copd+dys)
fit.0<-lmFit(gsva.val,mod0)
fit.0<-eBayes(fit.0)
res<-c()
res<-topTable(fit.0,coef=5:9,adjust.method="fdr",p.value=0.05,number=nrow(gsva.val))
write.table(res,file="dsyfact_wsex_smoke_copd_gsvares_fdr05.txt",sep="\t")
res<-topTable(fit.0,coef=5:9,adjust.method="fdr",p.value=1,number=nrow(gsva.val))
write.table(res,file="dsyfact_wsex_smoke_copd_gsvares.txt",sep="\t")

pdf(file="dys_wsmoke_sex_copd_fdr05_gsvavals.pdf")
sig.genes<-res$ID[which(res$adj.P.Val<0.05)]
data.heatmap<-exprs(gsva.val)[match(sig.genes,rownames(exprs(gsva.val))),]
heatmap3(data.heatmap,ColSideColors=as.matrix(t(column.matrix)),cexCol=0.2,cexRow=0.2,col.clustering = "unsupervised",col=colorRampPalette(c("blue","white","red"), space="rgb")(255))
dev.off()

dys.cont2<-as.character(annot.base$DYSGR_CAT)
dys.cont2<-gsub("Normal",0,dys.cont2)
dys.cont2<-gsub("Hyperplasia",0,dys.cont2)
dys.cont2<-gsub("Metaplasia",0,dys.cont2)
dys.cont2<-gsub("MildD",1,dys.cont2)
dys.cont2<-gsub("ModD",1,dys.cont2)
dys.cont2<-gsub("SevD",1,dys.cont2)
dys.cont2<-as.numeric(dys.cont2)
dys2<-as.factor(dys.cont2)
mod0<-model.matrix(~sex+smoke+copd+dys2)
fit.0<-lmFit(gsva.val,mod0)
fit.0<-eBayes(fit.0)
res<-c()
res<-topTable(fit.0,coef=5,adjust.method="fdr",p.value=1,number=nrow(gsva.val))
write.table(res,file="dsyfact_l2_wsex_smoke_copd_gsvares.txt",sep="\t")
pdf(file="dys_l2_wsmoke_sex_copd_fdr01_gsvavals_semi.pdf")
sig.genes<-res$ID[which(res$adj.P.Val<0.01)]
data.heatmap<-exprs(gsva.val)[match(sig.genes,rownames(exprs(gsva.val))),]
heatmap3(data.heatmap,ColSideColors=as.matrix(t(column.matrix)),cexCol=0.2,cexRow=0.2,col.clustering = "semisupervised",col=colorRampPalette(c("blue","white","red"), space="rgb")(255))
dev.off()
pdf(file="dys_l2_wsmoke_sex_copd_fdr001_gsvavals_semi.pdf")
sig.genes<-res$ID[which(res$adj.P.Val<0.001)]
data.heatmap<-exprs(gsva.val)[match(sig.genes,rownames(exprs(gsva.val))),]
heatmap3(data.heatmap,ColSideColors=as.matrix(t(column.matrix)),cexCol=0.2,cexRow=0.2,col.clustering = "semisupervised",col=colorRampPalette(c("blue","white","red"), space="rgb")(255))
dev.off()

#GSVA - try 2
data.comp.exp<-c();
data.comp.exp<-cpm.data.dge
entrez.id<-entrez[match(rownames(data.comp.exp),entrez[,1]),2]
d.entrez<-data.comp.exp[which(is.na(entrez.id)==FALSE),]
entrez.id<-entrez.id[which(is.na(entrez.id)==FALSE)]
red.id<-names(which(table(entrez.id[which(is.na(entrez.id)==FALSE)])>1))
nonred.id<-setdiff(entrez.id,red.id)
ind.final<-c()
for(i in 1:length(red.id)){
        ind<-which(entrez.id==red.id[i])
        mean.ind<-apply(d.entrez[ind,],1,mean)
        good.ind<-which(mean.ind==max(mean.ind))
        ind.final<-c(ind.final,ind[good.ind])
}
all.ind<-c(match(nonred.id,entrez.id),ind.final)
lab<-c(nonred.id,entrez.id[ind.final])
d.final<-d.entrez[all.ind,]
rownames(d.final)<-lab

pheno<-cbind(sex,copd,smoke)
rownames(pheno)<-colnames(d.final)
pheno.var<-c("sex","copd","smoke")
pheno.desc<-c("sex","copd","smoke")
pheno.df<-data.frame(pheno.var,pheno.desc)
pd <- new("AnnotatedDataFrame", data=data.frame(pheno),varMetadata=pheno.df)

e.set<-new("ExpressionSet",exprs=as.matrix(d.final),phenoData=pd,annotation="org.Hs.eg.db")
gsva.val <- gsva(e.set,gene.set,min.sz=10, max.sz=500,verbose=TRUE,rnaseq=TRUE)$es.obs
write.table(gsva.val,file="htseq_wcpmvals_GSVAvals.txt",sep="\t")

#Final model for now
mod0<-model.matrix(~sex+smoke+copd+dys)
fit.0<-lmFit(gsva.val,mod0)
fit.0<-eBayes(fit.0)
res<-c()
res<-topTable(fit.0,coef=5:9,adjust.method="fdr",p.value=0.05,number=nrow(gsva.val))
write.table(res,file="dsyfact_wsex_smoke_copd_gsvares_wcpmvals_fdr05.txt",sep="\t")
res<-topTable(fit.0,coef=5:9,adjust.method="fdr",p.value=1,number=nrow(gsva.val))
write.table(res,file="dsyfact_wsex_smoke_copd_gsvares_wcpmvals.txt",sep="\t")

pdf(file="dys_wsmoke_sex_copd_fdr05_gsvavals_wcpmvals.pdf")
sig.genes<-res$ID[which(res$adj.P.Val<0.05)]
data.heatmap<-exprs(gsva.val)[match(sig.genes,rownames(exprs(gsva.val))),]
heatmap3(data.heatmap,ColSideColors=as.matrix(t(column.matrix)),cexCol=0.2,cexRow=0.2,col.clustering = "unsupervised",col=colorRampPalette(c("blue","white","red"), space="rgb")(255))
dev.off()

#Get results using ENTREZ IDs
#create DGEList object
data.dge.ent<-DGEList(counts=d.final,genes=rownames(d.final))
#Scale Normalization
mod0<-model.matrix(~sex+smoke+copd+dys)
data.dge.ent<-calcNormFactors(data.dge.ent)
data.0<-voom(data.dge.ent,mod0,plot=FALSE)
fit.0<-lmFit(data.0,mod0)
fit.0<-eBayes(fit.0)
res<-c()
res<-topTable(fit.0,coef=5:9,adjust.method="fdr",p.value=1,number=nrow(d.final))
write.table(res,file="dysfact_wsmoke_sex_copd_ENTREZID.txt",sep="\t")

dys.cont2<-as.character(annot.base$DYSGR_CAT)
dys.cont2<-gsub("Normal",0,dys.cont2)
dys.cont2<-gsub("Hyperplasia",0,dys.cont2)
dys.cont2<-gsub("Metaplasia",1,dys.cont2)
dys.cont2<-gsub("MildD",2,dys.cont2)
dys.cont2<-gsub("ModD",2,dys.cont2)
dys.cont2<-gsub("SevD",2,dys.cont2)
dys.cont2<-as.numeric(dys.cont2)
dys2<-as.factor(dys.cont2)
mod0<-model.matrix(~sex+smoke+copd+dys2)
data.dge.ent<-calcNormFactors(data.dge.ent)
data.0<-voom(data.dge.ent,mod0,plot=FALSE)
fit.0<-lmFit(data.0,mod0)
fit.0<-eBayes(fit.0)
res<-c()
res<-topTable(fit.0,coef=5:6,adjust.method="fdr",p.value=1,number=nrow(d.final))
write.table(res,file="dysfact_wsmoke_sex_copd_l3_ENTREZID.txt",sep="\t")


dys.cont2<-as.character(annot.base$DYSGR_CAT)
dys.cont2<-gsub("Normal",0,dys.cont2)
dys.cont2<-gsub("Hyperplasia",0,dys.cont2)
dys.cont2<-gsub("Metaplasia",0,dys.cont2)
dys.cont2<-gsub("MildD",1,dys.cont2)
dys.cont2<-gsub("ModD",1,dys.cont2)
dys.cont2<-gsub("SevD",1,dys.cont2)
dys.cont2<-as.numeric(dys.cont2)
dys2<-as.factor(dys.cont2)
mod0<-model.matrix(~sex+smoke+copd+dys2)
data.dge.ent<-calcNormFactors(data.dge.ent)
data.0<-voom(data.dge.ent,mod0,plot=FALSE)
fit.0<-lmFit(data.0,mod0)
fit.0<-eBayes(fit.0)
res<-c()
res<-topTable(fit.0,coef=5,adjust.method="fdr",p.value=1,number=nrow(d.final))
write.table(res,file="dysfact_wsmoke_sex_copd_l2_ENTREZID.txt",sep="\t")

