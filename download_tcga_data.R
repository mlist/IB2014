#All data used in this publication can be downloaded from
#https://tcga-data.nci.nih.gov/docs/publications/brca_2012/
  
#download gene expression data with subtype information
brca.rna.dest <- "data/BRCA.exp.547.med.txt"
brca.pam50 <- "data/BRCA.547.PAM50.SigClust.Subtypes.txt"
download.file("http://tcga-data.nci.nih.gov/docs/publications/brca_2012/BRCA.exp.547.med.txt", destfile=brca.rna.dest, method="wget")
download.file("http://tcga-data.nci.nih.gov/docs/publications/brca_2012/BRCA.547.PAM50.SigClust.Subtypes.txt", destfile=brca.pam50, method="wget")

#download methylation data
temp <- tempfile()
download.file("http://tcga-data.nci.nih.gov/docs/publications/brca_2012/BRCA.methylation.27k.450k.zip", destfile=temp, method="wget")
unzip(temp, "BRCA.methylation.27k.450k.txt", exdir="data")
unlink(temp)

#download gene mapping data for methylation data
temp <- tempfile()
download.file("https://tcga-data.nci.nih.gov/docs/publications/tcga/integration/adfs/tcga/jhu-usc.edu_TCGA_HumanMethylation27.v2.adf.txt.zip", destfile=temp, method="wget")
unzip(temp, "jhu-usc.edu_TCGA_HumanMethylation27.v2.adf.txt", exdir="data")
unlink(temp)

