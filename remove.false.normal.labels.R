#here we take the output of process_tcga_data.R and remove the 8 normal samples that PAM50 mis-labelled
which.samples <- names(PAM50)[which(PAM50 != "Normal")]
PAM50 <- as.factor(as.character(PAM50[which.samples]))
names(PAM50) <- which.samples
SigClust <- SigClust[which.samples]
names(SigClust) <- which.samples

brca.gene.expr <- brca.gene.expr[which.samples,]
brca.meth.full <- brca.meth.full[which.samples,]
combinedFeatureMatrix <- combinedFeatureMatrix[which.samples,]

save(list=c("brca.gene.expr", "brca.meth.full", "combinedFeatureMatrix", "brca.subtypes", "PAM50", "SigClust"), file="data/tcga_gene_expr_and_meth_no_normal.RData")

#PAM50 random forest no normal 
library(limma)
PAM50.genes <- read.delim("../lists//pam50.txt")
PAM50.genes <- sapply(PAM50.genes$x, as.character)
PAM50.rf <- randomForest(brca.gene.expr[,PAM50.genes], PAM50, ntree=5000, mtry=10)
save(PAM50, PAM50.rf, file="../results_no_normal/PAM50.randomForest.RData")
