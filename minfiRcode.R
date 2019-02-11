source("https://bioconductor.org/biocLite.R")
biocLite("minfi")
library('minfi')
#Read Idat files
idatPath <- ("C:/Users/swapnajoshi-admin/Documents/GIT_workspace/DNA_meth/DNA_methylation_02062017/Methylumi_analysis/Idat_all")

list.files(idatPath)
target1 <- read.delim("C:/Users/swapnajoshi-admin/Documents/GIT_workspace/DNA_meth/DNA_methylation_02062017/Methylumi_analysis/Idat_all/SampleSheet.txt", as.is = TRUE); dim(target1)
#names(target1)
#target1$Basename
target1$Basename <- file.path(idatPath, target1$Basename)
rgset <- read.450k(target1$Basename, verbose = TRUE)

row.names(target1) <- substr(target1$Basename)
pData(rgset) <- target1
save(rgset, file = "C:/Users/swapnajoshi-admin/Documents/GIT_workspace/DNA_meth/DNA_methylation_02062017/Minfi_analysis/rgset.rda")
dim(getRed(rgset))
#[1] 622399    288
dim(getGreen(rgset))
#[1] 622399    288

source("https://bioconductor.org/biocLite.R")
biocLite("IlluminaHumanMethylation450kmanifest")
source("https://bioconductor.org/biocLite.R")
biocLite("IlluminaHumanMethylation450kanno.ilmn12.hg19")

rgset_proc <- preprocessFunnorm(rgset)
pCutoff <- 0.01
pvals <- detectionP(rgset)
is.na(assays(rgset_proc)$Beta) <- (pvals[rownames(rgset_proc), colnames(rgset_proc)] >= pCutoff)

save(rgset_proc, file = "C:/Users/swapnajoshi-admin/Documents/GIT_workspace/DNA_meth/DNA_methylation_02062017/Minfi_analysis/rgset_proc.rda")

granges(rgset_proc)

getBeta(rgset_proc)[1:3,1:3]


rgset_beta <- getBeta(rgset_proc)
save(rgset_beta, file = "C:/Users/swapnajoshi-admin/Documents/GIT_workspace/DNA_meth/DNA_methylation_02062017/Minfi_analysis/rgset_beta.rda")

rgset_proc_pbmc <- rgset_proc[,grep("_P",pData(rgset_proc)$External.Sample.ID)]
rgset_proc_col <- rgset_proc[,grep("_B",pData(rgset_proc)$External.Sample.ID)]
pbmc_beta <- getBeta(rgset_proc_pbmc)
col_beta <- getBeta(rgset_proc_col)
save(pbmc_beta, file = "C:/Users/swapnajoshi-admin/Documents/GIT_workspace/DNA_meth/DNA_methylation_02062017/Minfi_analysis/pbmc_beta.rda")
save(col_beta, file = "C:/Users/swapnajoshi-admin/Documents/GIT_workspace/DNA_meth/DNA_methylation_02062017/Minfi_analysis/col_beta.rda")

load("C:/Users/swapnajoshi-admin/Documents/GIT_workspace/DNA_meth/DNA_methylation_02062017/Minfi_analysis/rgset.rda")
library(shinyMethyl)
summaryP <- shinySummarize(rgset)
runShinyMethyl(summaryP)



load("C:/Users/swapnajoshi-admin/Documents/GIT_workspace/DNA_meth/DNA_methylation_02062017/Minfi_analysis/rgset.rda")
rgset_pbmc <- rgset[,grep("_P",pData(rgset)$External.Sample.ID)]
rgset_col <- rgset[,grep("_B",pData(rgset)$External.Sample.ID)]
GRset.quantile.pbmc <- preprocessQuantile(rgset_pbmc, fixOutliers = TRUE,
  removeBadSamples = TRUE, badSampleCutoff = 10.5,
  quantileNormalize = TRUE, stratified = TRUE, 
  mergeManifest = FALSE, sex = NULL)
save(GRset.quantile.pbmc, file = "C:/Users/swapnajoshi-admin/Documents/GIT_workspace/DNA_meth/DNA_methylation_02062017/Minfi_analysis/GRset.quantile.pbmc.rda")
GRset.quantile.col <- preprocessQuantile(rgset_col, fixOutliers = TRUE,
  removeBadSamples = TRUE, badSampleCutoff = 10.5,
  quantileNormalize = TRUE, stratified = TRUE, 
  mergeManifest = FALSE, sex = NULL)
save(GRset.quantile.col, file = "C:/Users/swapnajoshi-admin/Documents/GIT_workspace/DNA_meth/DNA_methylation_02062017/Minfi_analysis/GRset.quantile.col.rda")
library(FlowSorted.Blood.450k)
cellCounts <- estimateCellCounts(RGset.pbmc)
