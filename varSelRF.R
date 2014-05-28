library(snow)
library(varSelRF)

#load data
load("~/BCsubtypes/data/tcga_gene_expr_and_meth.RData")

#start cluster, executing in parallel
cl <- makeMPIcluster(16) #uses 16 cores. change to number of available cores.

#perform normal varSelRF runs for each model
brca.gene.expr.varSelRF <- varSelRF(brca.gene.expr, PAM50, keep.forest=T)
save(brca.gene.expr.varSelRF, file="~/BCsubtypes/results/brca.gene.expr.varSelRF.RData")

brca.meth.full.varSelRF <- varSelRF(brca.meth.full, PAM50, keep.forest=T)
save(brca.meth.full.varSelRF, file="~/BCsubtypes/results/brca.meth.full.varSelRF.RData")

brca.combined.selrf <- varSelRF(combinedFeatureMatrix, PAM50, keep.forest=T)
save(brca.combined.selrf, file="~/BCsubtypes/results/brca.combined.selrf.RData")

#bootstrapping on each model using the above results as seed.
brca.gene.expr.boot <- varSelRFBoot(brca.gene.expr, PAM50, bootnumber=10, usingCluster=T, TheCluster=cl, srf=brca.gene.expr.varSelRF)
save(brca.gene.expr.boot, file="~/BCsubtypes/results/brca.gene.expr.selrf.boot.RData")

brca.meth.full.boot <- varSelRFBoot(brca.meth.full, PAM50, bootnumber=10, usingCluster=T, TheCluster=cl, srf=brca.meth.full.varSelRF)
save(brca.meth.full.boot, file="~/BCsubtypes/results/brca.meth.full.selrf.boot.RData")

brca.combined.selrf.boot <- varSelRFBoot(combinedFeatureMatrix, PAM50, bootnumber=10, usingCluster=T, TheCluster=cl, srf=brca.combined.selrf)
save(brca.combined.selrf.boot, file="~/BCsubtypes/results/brca.combined.selrf.boot.RData")

#stop cluster
stopCluster(cl)
