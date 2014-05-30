These are the R scripts used in 

List et al. Classification of Breast Cancer Subtypes by Combining Gene Expression and DNA Methylation Data, Journal of Integrative Bioinformatics (in press)

This work was presented at the 10th international symposium on Integrative Bioinformatics on 12th of May 2014 in Newcastle, UK.

The following files are included:
* download_tcga_data.R: script for downloading the necessary files from the TCGA website
* process_tcga_data.R: script to process the downloaded TCGA breast cancer data sets
* remove.false.normal.labels.R: script to remove samples designated as normal by PAM50 from the processed data. helpful if they should be excluded from the analysis.
* merge.normal.samples.R: Add non-tumor "normal" samples 
* varSelRF.R: Perform random forest using the varSelRF package.
* helper_methods.R: Neat methods for analysing the data
* pam50.txt: The 50 genes that are part of the PAM50 classifier (gold standard).
