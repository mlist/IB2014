#### Gene Expression Data and subtype labels ####

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

#### Methylation Data ####

#load full methylation data set
brca.meth.full.raw <- read.delim("data/BRCA.methylation.27k.450k.txt")

#probe to gene symbol
meth.probe.map <- read.delim("data/jhu-usc.edu_TCGA_HumanMethylation27.v2.adf.txt")
meth.probe.id.to.gene.symbol <- meth.probe.map$SYMBOL
names(meth.probe.id.to.gene.symbol) <- meth.probe.map$IlmnID

#how many samples can we match to the ones we know the subtype for? filter step 2
commonSamples <- intersect(colnames(brca.meth.full.raw),substr(names(PAM50),1,16))
cat(paste("Number of samples we have subtype information on:", length(commonSamples) , "/", ncol(brca.meth.full.raw), "\n"))
brca.meth.full <- brca.meth.full.raw[,commonSamples]

#replace probe names with gene symbols (better do this on the feature set, takes not so much time and we can distuinguish duplicates)
meth.probe.genes <- meth.probe.id.to.gene.symbol[row.names(brca.meth.full)]
unique.meth.probe.genes <- levels(meth.probe.genes)[-1] #omitting first entry of unique.meth.probe.genes because its the empty string

# take median of gene duplicates
# extract rows for all genes and merge them into one, then bind rows together
# NOTE: when running the combined data of gene expression and methylation data the features of the methylation data 
# needs to be marked to separate the gene from the meth features in the results. Alternatively you can run random forest
# with probe names (aggregated) as feature names to make sure you recognize meth features. However, since the selected features from 
# the combined matrix exclusively consisted of meth features in all runs, we use gene identifiers here for publication.
library(foreach)
brca.meth.full.merged <- foreach(unique.gene=unique.meth.probe.genes, .combine=rbind) %do% 
{
  probes.for.this.gene <- names(meth.probe.genes[meth.probe.genes == unique.gene])
  apply(brca.meth.full[probes.for.this.gene,],2, median, na.rm=T)
}
row.names(brca.meth.full.merged) <- unique.meth.probe.genes

brca.meth.full <- brca.meth.full.merged
rm(brca.meth.full.merged)

#remove features that are not complete, e.g. contain NAs. These are not allowed in random forest. filter step 3
cat(paste("Number of probes after omitting probes with NA values (out of", nrow(brca.meth.full.raw), "):"))
brca.meth.full <- brca.meth.full[complete.cases(brca.meth.full),]
cat(paste(nrow(brca.meth.full), "\n"))

#transpose
brca.meth.full <- t(brca.meth.full)

#remove samples from gene expression data that is not also found in methylation data
row.names(brca.gene.expr) <- substr(row.names(brca.gene.expr),1,16)
brca.gene.expr <- brca.gene.expr[commonSamples,]
names(PAM50) <- substr(names(PAM50),1,16)
PAM50 <- PAM50[commonSamples]

names(SigClust) <- substr(names(SigClust),1,16)
SigClust <- SigClust[commonSamples]

#check that all sample names match
all.equal(row.names(brca.gene.expr), names(PAM50))
all.equal(row.names(brca.meth.full), names(PAM50))

combinedFeatureMatrix <- cbind(brca.gene.expr, brca.meth.full)

save(list=c("brca.gene.expr", "brca.meth.full", "combinedFeatureMatrix", "brca.subtypes", "PAM50", "SigClust"), file="data/tcga_gene_expr_and_meth.RData")


