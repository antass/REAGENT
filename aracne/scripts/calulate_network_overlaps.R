library(org.Hs.eg.db)

# convert easily from affy probe ID to EntrezGene ID
x <- org.Hs.egSYMBOL
# Get the gene symbol that are mapped to an entrez gene identifiers
mapped_genes <- mappedkeys(x)
# Convert to a list
xx <- as.list(x[mapped_genes])

# load networks
load("~/Meta_Analysis/kyoto/aracne/data/airway_FS2_ARACNE.rda")  # air.aracne.fs2.network
load("~/Meta_Analysis/kyoto/aracne/data/allegro_FS2_ARACNE.rda")  # all.aracne.fs2.network
load("~/Meta_Analysis/kyoto/aracne/data/airway_cancer_FS2_ARACNE.rda")  # air.cancer.aracne.fs2.network
load("~/Meta_Analysis/kyoto/aracne/data/airway_normal_FS2_ARACNE.rda")  # air.normal.aracne.fs2.network
load("~/Meta_Analysis/kyoto/aracne/data/allegro_cancer_FS2_ARACNE.rda")  # all.cancer.aracne.fs2.network
load("~/Meta_Analysis/kyoto/aracne/data/allegro_normal_FS2_ARACNE.rda")  # all.normal.aracne.fs2.network
lnames<-load("~/Meta_Analysis/kyoto/aracne/data/allegro_ARACNE.rda")  # all.network
lnames<-load("~/Meta_Analysis/kyoto/aracne/data/airway_ARACNE.rda")  # air.network






# # convert Entrez Gene IDs to Gene Symbols
# entrez <- sapply(colnames(air.aracne.fs2.network), function(i) substr(i, 1, nchar(i)-3))
# symbols <- xx[entrez]
# dimnames(air.aracne.fs2.network) <- list(symbols, symbols)
# dimnames(all.aracne.fs2.network) <- list(symbols, symbols)
# dimnames(air.cancer.aracne.fs2.network) <- list(symbols, symbols)
# dimnames(air.normal.aracne.fs2.network) <- list(symbols, symbols)
# dimnames(all.cancer.aracne.fs2.network) <- list(symbols, symbols)
# dimnames(all.normal.aracne.fs2.network) <- list(symbols, symbols)
# dimnames(air.network) <- list(symbols, symbols)
# dimnames(all.network) <- list(symbols, symbols)
# 
# save(air.aracne.fs2.network, file="~/Meta_Analysis/kyoto/aracne/data/airway_symbol_FS2_ARACNE.rda")  # air.aracne.fs2.network
# save(all.aracne.fs2.network, file="~/Meta_Analysis/kyoto/aracne/data/allegro_symbol_FS2_ARACNE.rda")  # all.aracne.fs2.network
# save(air.cancer.aracne.fs2.network, file="~/Meta_Analysis/kyoto/aracne/data/airway_cancer_symbol_FS2_ARACNE.rda")  # air.cancer.aracne.fs2.network
# save(air.normal.aracne.fs2.network, file="~/Meta_Analysis/kyoto/aracne/data/airway_normal_symbol_FS2_ARACNE.rda")  # air.normal.aracne.fs2.network
# save(all.cancer.aracne.fs2.network, file="~/Meta_Analysis/kyoto/aracne/data/allegro_cancer_symbol_FS2_ARACNE.rda")  # all.cancer.aracne.fs2.network
# save(all.normal.aracne.fs2.network, file="~/Meta_Analysis/kyoto/aracne/data/allegro_normal_symbol_FS2_ARACNE.rda")  # all.normal.aracne.fs2.network
# save(air.network, file="~/Meta_Analysis/kyoto/aracne/data/airway_symbol_ARACNE.rda")  # air.aracne.fs2.network
# save(all.network, file="~/Meta_Analysis/kyoto/aracne/data/allegro_symbol_ARACNE.rda")  # all.aracne.fs2.network



# convert weights to binary (edge exists or not)
air.ar.vector <- as.vector(air.aracne.fs2.network)
air.ar.vector <- ifelse(air.ar.vector==0, 0, 1)  # replace weights with binary (if mi > 0 then mi = 1)
air.ar.matrix <- matrix(air.ar.vector, nrow=nrow(air.aracne.fs2.network), byrow=FALSE, dimnames=dimnames(air.aracne.fs2.network))

air.c.ar.vector <- as.vector(air.cancer.aracne.fs2.network)
air.c.ar.vector <- ifelse(air.c.ar.vector==0, 0, 1)  # replace weights with binary (if mi > 0 then mi = 1)
air.c.ar.matrix <- matrix(air.c.ar.vector, nrow=nrow(air.cancer.aracne.fs2.network), byrow=FALSE, dimnames=dimnames(air.cancer.aracne.fs2.network))

air.n.ar.vector <- as.vector(air.normal.aracne.fs2.network)
air.n.ar.vector <- ifelse(air.n.ar.vector==0, 0, 1)  # replace weights with binary (if mi > 0 then mi = 1)
air.n.ar.matrix <- matrix(air.n.ar.vector, nrow=nrow(air.normal.aracne.fs2.network), byrow=FALSE, dimnames=dimnames(air.normal.aracne.fs2.network))

air.orig.ar.vector <- as.vector(air.aracne)
air.orig.ar.vector <- ifelse(air.orig.ar.vector==0, 0, 1)  # replace weights with binary (if mi > 0 then mi = 1)
air.orig.ar.matrix <- matrix(air.orig.ar.vector, nrow=nrow(air.aracne), byrow=FALSE, dimnames=dimnames(air.aracne))


all.ar.vector <- as.vector(all.aracne.fs2.network)
all.ar.vector <- ifelse(all.ar.vector==0, 0, 1)  # replace weights with binary (if mi > 0 then mi = 1)
all.ar.matrix <- matrix(all.ar.vector, nrow=nrow(all.aracne.fs2.network), byrow=FALSE, dimnames=dimnames(all.aracne.fs2.network))

all.c.ar.vector <- as.vector(all.cancer.aracne.fs2.network)
all.c.ar.vector <- ifelse(all.c.ar.vector==0, 0, 1)  # replace weights with binary (if mi > 0 then mi = 1)
all.c.ar.matrix <- matrix(all.c.ar.vector, nrow=nrow(all.cancer.aracne.fs2.network), byrow=FALSE, dimnames=dimnames(all.cancer.aracne.fs2.network))

all.n.ar.vector <- as.vector(all.normal.aracne.fs2.network)
all.n.ar.vector <- ifelse(all.n.ar.vector==0, 0, 1)  # replace weights with binary (if mi > 0 then mi = 1)
all.n.ar.matrix <- matrix(all.n.ar.vector, nrow=nrow(all.normal.aracne.fs2.network), byrow=FALSE, dimnames=dimnames(all.normal.aracne.fs2.network))

all.orig.ar.vector <- as.vector(all.aracne)
all.orig.ar.vector <- ifelse(all.orig.ar.vector==0, 0, 1)  # replace weights with binary (if mi > 0 then mi = 1)
all.orig.ar.matrix <- matrix(all.orig.ar.vector, nrow=nrow(all.aracne), byrow=FALSE, dimnames=dimnames(all.aracne))


# identify hubs
cons <- apply(air.ar.matrix, 2, sum)
cons.air.full <- cons[order(cons, decreasing=TRUE)] 

cons <- apply(air.c.ar.matrix, 2, sum)
cons.air.cancer <- cons[order(cons, decreasing=TRUE)] 

cons <- apply(air.n.ar.matrix, 2, sum)
cons.air.normal <- cons[order(cons, decreasing=TRUE)] 

cons <- apply(air.orig.ar.matrix, 2, sum)
cons.air.original <- cons[order(cons, decreasing=TRUE)] 

cons <- apply(all.ar.matrix, 2, sum)
cons.all.full <- cons[order(cons, decreasing=TRUE)] 

cons <- apply(all.c.ar.matrix, 2, sum)
cons.all.cancer <- cons[order(cons, decreasing=TRUE)] 

cons <- apply(all.n.ar.matrix, 2, sum)
cons.all.normal <- cons[order(cons, decreasing=TRUE)] 

cons <- apply(all.orig.ar.matrix, 2, sum)
cons.all.original <- cons[order(cons, decreasing=TRUE)] 

neighbors <- function(id, network){
  cons <- network[rownames(network)==id, ]
  pos <- cons[cons==1]
  genes <- names(pos)
  return(genes)
}

nbr.air.c <- neighbors(names(cons.air.cancer[1]), air.c.ar.matrix)
nbr.air.n <- neighbors(names(cons.air.normal[1]), air.n.ar.matrix)
nbr.air.f <- neighbors(names(cons.air.full[1]), air.ar.matrix)
nbr.air.o <- neighbors(names(cons.air.original[1]), air.orig.ar.matrix)

nbr.all.c <- neighbors(names(cons.all.cancer[1]), all.c.ar.matrix)
nbr.all.n <- neighbors(names(cons.all.normal[1]), all.n.ar.matrix)
nbr.all.f <- neighbors(names(cons.all.full[1]), all.ar.matrix)
nbr.all.o <- neighbors(names(cons.all.original[1]), all.orig.ar.matrix)

genes <- rownames(air.c.ar.matrix)

df <- data.frame()
for(gene in genes){
  c <- neighbors(gene, air.c.ar.matrix)
  n <- neighbors(gene, air.n.ar.matrix)
  f <- neighbors(gene, air.ar.matrix)
  x <- data.frame(gene=gene, cancer=length(c), normal=length(n), fs=length(f), cn=length(intersect(c,n)), 
                  cf=length(intersect(c,f)), nf=length(intersect(n,f)))
#   print(x)
  df <- rbind(df, x)
}

df[order(df$cancer, decreasing=TRUE),]

### Make figure like in paper : # of interactions per particular strength
air.c.ar.vector.interactions <- as.vector(air.cancer.aracne.fs2.network)
table(air.c.ar.vector.interactions)
summary(air.c.ar.vector.interactions[air.c.ar.vector.interactions!=0])
xc <- air.c.ar.vector.interactions[air.c.ar.vector.interactions!=0]

air.n.ar.vector.interactions <- as.vector(air.normal.aracne.fs2.network)
table(air.n.ar.vector.interactions)
summary(air.n.ar.vector.interactions[air.n.ar.vector.interactions!=0])
xn <- air.n.ar.vector.interactions[air.n.ar.vector.interactions!=0]


xxc <- cut(xc, breaks=seq(0,max(xc), by=0.001), include.lowest=FALSE)
table(xxc)
xcdf <- data.frame(cbind(bins=names(table(xxc)), interactions=as.numeric(as.vector(table(xxc))), weight=as.numeric(seq(0,max(xc), by=0.001)[1:length(seq(0, max(xc), by=0.001))-1])))


xxn <- cut(xn, breaks=seq(0,max(xn), by=0.001), include.lowest=FALSE)
table(xxn)
xndf <- data.frame(cbind(bins=names(table(xxn)), interactions=as.numeric(as.vector(table(xxn))), weight=as.numeric(seq(0,max(xn), by=0.001)[1:length(seq(0, max(xn), by=0.001))-1])))




### remove redundant edges
x <- air.cancer.aracne.fs2.network
upperTriangle(x, diag=FALSE) <- NA
xx <- melt(x)
xx <- xx[!is.na(xx[,3]) & !(xx[,3]==0),]
wgts <- xx[,3]
bins <- cut(wgts, seq(0, max(wgts), by=0.01), include.lowest=TRUE)
xcdf <- data.frame(cbind(bins=names(table(bins)), interactions=as.numeric(as.vector(table(bins))), weight=as.numeric(seq(0, max(wgts), by=0.01)[1:length(seq(0, max(wgts), by=0.01))-1])))
xcdf <- xcdf[ xcdf$interactions!=0,]

x <- air.normal.aracne.fs2.network
upperTriangle(x, diag=FALSE) <- NA
xx <- melt(x)
xx <- xx[!is.na(xx[,3]) & !(xx[,3]==0),]
wgts <- xx[,3]
bins <- cut(wgts, seq(0, max(wgts), by=0.01), include.lowest=TRUE)
xndf <- data.frame(cbind(bins=names(table(bins)), interactions=as.numeric(as.vector(table(bins))), weight=as.numeric(seq(0, max(wgts), by=0.01)[1:length(seq(0, max(wgts), by=0.01))-1])))
xndf <- xcdf[ xndf$interactions!=0,]

x <- all.cancer.aracne.fs2.network
upperTriangle(x, diag=FALSE) <- NA
xx <- melt(x)
xx <- xx[!is.na(xx[,3]) & !(xx[,3]==0),]
wgts <- xx[,3]
bins <- cut(wgts, seq(0, max(wgts), by=0.01), include.lowest=TRUE)
xcdf <- data.frame(cbind(bins=names(table(bins)), interactions=as.numeric(as.vector(table(bins))), weight=as.numeric(seq(0, max(wgts), by=0.01)[1:length(seq(0, max(wgts), by=0.01))-1])))
xcdf <- xcdf[ xcdf$interactions!=0,]

x <- all.normal.aracne.fs2.network
upperTriangle(x, diag=FALSE) <- NA
xx <- melt(x)
xx <- xx[!is.na(xx[,3]) & !(xx[,3]==0),]
wgts <- xx[,3]
bins <- cut(wgts, seq(0, max(wgts), by=0.01), include.lowest=TRUE)
xndf <- data.frame(cbind(bins=names(table(bins)), interactions=as.numeric(as.vector(table(bins))), weight=as.numeric(seq(0, max(wgts), by=0.01)[1:length(seq(0, max(wgts), by=0.01))-1])))
xndf <- xcdf[ xndf$interactions!=0,]


## MAKE DIAGNOSTIC PLOTS
# interaction N vs strength
# fit a smooth curve

# airway or allegro (adjust above)
ints.c <- as.numeric(as.character(xcdf$interactions))
ints2.c <- ints.c[ints.c!=0]
wgts.c <- as.numeric(as.character(xcdf$weight))[ints.c!=0]

lo.c <- loess(wgts.c ~ ints2.c)
plot(x=ints2.c, y=wgts.c)
j.c <- order(ints2.c)
lines(x=ints2.c[j.c], lo.c$fitted[j.c], col="red", lwd=3)

ints.n <- as.numeric(as.character(xndf$interactions))
ints2.n <- ints.n[ints.n!=0]
wgts.n <- as.numeric(as.character(xndf$weight))[ints.n!=0]

lo.n <- loess(wgts.n ~ ints2.n)
plot(x=ints2.n, y=wgts.n)
j.n <- order(ints2.n)
lines(x=ints2.n[j.n], lo.n$fitted[j.n], col="red", lwd=3)

pdf("~/Meta_Analysis/kyoto/aracne/figures/Allegro_ARACNE_interactions_cancer_normal.pdf")
png("~/Meta_Analysis/kyoto/aracne/figures/Allegro_ARACNE_interactions_cancer_normal.png")
pdf("~/Meta_Analysis/kyoto/aracne/figures/Airway_ARACNE_interactions_cancer_normal.pdf")
png("~/Meta_Analysis/kyoto/aracne/figures/Airway_ARACNE_interactions_cancer_normal.png")
plot(x=ints2.c[j.c], lo.c$fitted[j.c], col="red", lwd=3, type="line", 
     ylab="Interaction weights", xlab="Number of interactions") 
title("Interaction strength in Allegro\n cancer and normal networks")
lines(x=ints2.n[j.n], lo.n$fitted[j.n], col="blue", lwd=3)
legend(30, 0.8, c("normal", "cancer"), col=c("blue", "red"), lwd=3)
dev.off()




