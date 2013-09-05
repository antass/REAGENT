# example data frame
set.seed(33)
expr <- data.frame(EnsemblID=c("A", "B", "B", "C", "C", "C", "A", "B", "A", "C", "D"),
                   sample1=rnorm(11),
                   sample2=rnorm(11),
                   sample3=rnorm(11))


probes2genes <- function(expr) {
  expr2 <- aggregate(expr[,-1], by=list(Gene=as.factor(expr$EnsemblID)), FUN=median)
  print(dim(expr))
  print(dim(expr2))
  
  good.index <- c()
  for (profile in seq_along(expr$EnsemblID)) {
    if(profile%%100==0) {print(profile)}
    gene <- as.numeric(expr[profile, 2:ncol(expr)])
    median <- as.numeric(expr2[expr2$Gene == expr$EnsemblID[profile], 2:ncol(expr2)])
#     print(gene)
#     print(median)
#     print(cor.test(gene, median)$estimate)
    if (cor.test(gene, median)$estimate > 0.5) {  # just >0.7 because I want the correlation to act one direction
      good.index <- c(good.index, profile)
    }
  }
  
  print("second round")
  expr <- expr[good.index, ]
  
  expr2 <- aggregate(expr[,-1], by=list(Gene=as.factor(expr$EnsemblID)), FUN=median)
  
  expr$EnsemblID
  good.index <- c()
  for (profile in seq_along(expr$EnsemblID)) {
    gene <- as.numeric(expr[profile, 2:ncol(expr)])
    median <- as.numeric(expr2[expr2$Gene == expr$EnsemblID[profile], 2:ncol(expr2)])
#     print(gene)
#     print(median)
#     print(cor.test(gene, median)$estimate)
    if (cor.test(gene, median)$estimate > 0.5) {  # just >0.7 because I want the correlation to act one direction
      good.index <- c(good.index, profile)
    }
  }
  
  print(nrow(expr)==length(good.index))
  return(expr2)
}

expr.final <- probes2genes(expr)