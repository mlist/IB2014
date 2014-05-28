#in this script we take the methylation data output from process_tcga_data.R and add the normal samples to the methylation matrix
#gene expression data is created from scratch

#load gene expression and subtype data
brca.gene.expr <- read.delim("data/BRCA.exp.547.med.txt", sep="\t", header=T, row.names=1)
brca.subtypes <- read.delim("data/BRCA.547.PAM50.SigClust.Subtypes.txt", header=T, row.names=1)

#filter for tumor samples and repair sample names 
library(stringr)
row.names(brca.subtypes) <- str_replace_all(row.names(brca.subtypes), "-", ".") 

#select which samples to filter for
#brca.subtypes <- subset(brca.subtypes, Type=="tumor")
brca.gene.expr <- brca.gene.expr[,colnames(brca.gene.expr) %in% row.names(brca.subtypes)]

#order labels according to sample ordering of data
brca.subtypes <- brca.subtypes[order(row.names(brca.subtypes), colnames(brca.gene.expr)),]

#check if ordering worked
if(!all(row.names(brca.subtypes) == colnames(brca.gene.expr))) stop("ordering of sample labels failed!!!!")

#omit genes with NA values
brca.gene.expr <- na.omit(brca.gene.expr)

#transform gene expression matrix such that the genes are the columns and the rows correspond to patients
brca.gene.expr <- t(brca.gene.expr)

#use PAM50 classifier as class labels
PAM50 <- brca.subtypes[,c("PAM50")]
names(PAM50) <- row.names(brca.subtypes)

#Alternativ SigClust
SigClust <- brca.subtypes[,c("Siglust")]
names(SigClust) <- row.names(brca.subtypes)

brca.meth.full.raw <- read.delim("data/BRCA.methylation.27k.450k.txt")

normal.samples <- intersect(colnames(brca.meth.full.raw), substr(row.names(subset(brca.subtypes, PAM50=="Normal" & Type!="tumor")), 1,16))
brca.meth.full.normal  <- brca.meth.full.raw[,normal.samples]
#replace probe names with gene symbols (better do this on the feature set, takes not so much time and we can distuinguish duplicates)
meth.probe.genes <- meth.probe.id.to.gene.symbol[row.names(brca.meth.full.normal)]
unique.meth.probe.genes <- levels(meth.probe.genes)[-1] #omitting first entry of unique.meth.probe.genes because its the empty string

# take median of gene duplicates
library(foreach)
brca.meth.full.normal.merged <- foreach(unique.gene=unique.meth.probe.genes, .combine=rbind) %do% 
{
  probes.for.this.gene <- names(meth.probe.genes[meth.probe.genes == unique.gene])
  apply(brca.meth.full.normal[probes.for.this.gene,],2, median, na.rm=T)
}
row.names(brca.meth.full.normal.merged) <- unique.meth.probe.genes

brca.meth.full.normal <- brca.meth.full.normal.merged
rm(brca.meth.full.normal.merged)

brca.meth.full.normal <- brca.meth.full.normal[complete.cases(brca.meth.full.normal),]
brca.meth.full.normal <- t(brca.meth.full.normal)
#find common columns 
commonCols <- intersect(colnames(brca.meth.full.normal), colnames(brca.meth.full))

brca.meth.full.normal.merged <- rbind(brca.meth.full.normal[, commonCols], brca.meth.full[, commonCols])

commonSamples <- intersect(rownames(brca.meth.full.normal.merged),substr(names(PAM50),1,16))

#restore correct sample order
brca.meth.full.normal <- brca.meth.full.normal.merged[commonSamples,]

#remove samples from gene expression data that is not also found in methylation data
row.names(brca.gene.expr) <- substr(row.names(brca.gene.expr),1,16)
brca.gene.expr <- brca.gene.expr[commonSamples,]
names(PAM50) <- substr(names(PAM50),1,16)
PAM50 <- PAM50[commonSamples]

names(SigClust) <- substr(names(SigClust),1,16)
SigClust <- SigClust[commonSamples]

#check that all sample names match
all.equal(row.names(brca.gene.expr), names(PAM50))
all.equal(row.names(brca.meth.full.normal), names(PAM50))

brca.meth.full <- brca.meth.full.normal 
combinedFeatureMatrix <- cbind(brca.gene.expr, brca.meth.full.normal)

save(list=c("brca.gene.expr", "brca.meth.full", "combinedFeatureMatrix", "PAM50", "SigClust"), file="data/tcga_gene_expr_and_meth_incl_normal.RData")
