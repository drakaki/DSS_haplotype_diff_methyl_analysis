library(DSS)
require(bsseq)
path_dss = file.path(system.file(package="DSS"), "extdata")
path="/path/to/bed/files"
setwd(path)

# get hap1 replicates
ext1=".cpg.hap1.reference.mincov4.bed"
hap1.files = list.files(path=path, pattern=ext1)

counter=0
for (file in hap1.files){
    counter=counter+1
    file_name <- paste0("hap1.",counter)
    assign(file_name,read_table(file,col_names = FALSE))
  
}

# get hap2 replicatess
ext2=".cpg.hap2.reference.mincov4.bed"
hap2.files = list.files(path=path, pattern=ext2)

counter=0
for (file in hap2.files){
  counter=counter+1
  file_name <- paste0("hap2.",counter)
  assign(file_name,read_table(file,col_names = FALSE))
  
}

hap_list <- mget(grep("hap[1-2].[1-6]", ls(), value = TRUE))
# filtering
hap_list <- lapply(hap_list, function(x) x %>% dplyr::select(X1,X2,X6,X7) %>% filter(!(X1=="chrX")))
# colnames as instructed in DSS package vignette
cols <- c("chr","pos","N","X")
hap_list <- lapply(hap_list, setNames, cols)


# Replicate 3 has too few CpG sites covered
which.min(unlist(lapply(hap_list, function(x) nrow(x)))) #hap2.3
# split BSobj to save memory usage
Bsobj_1 <- makeBSseqData(hap_list[grep("hap1",names(hap_list))[-3]], 
                         grep("hap1",names(hap_list), value=T)[-3]) #sampleNames
BSobj_2 <- makeBSseqData(hap_list[grep("hap2",names(hap_list))[-3]], 
                         grep("hap2",names(hap_list), value=T)[-3]) #sampleNames

BSobj <- combine(BSobj_1, BSobj_2)

cov <- as.data.frame(getCoverage(BSobj))
cov$pos  <- BSobj@rowRanges@ranges@start
k1 <- which(rowSums(cov[,grep("hap1", colnames(cov)])>=4)
k2 <- which(rowSums(cov[,grep("hap2",colnames(cov))])>=4)
keep_sites <- Reduce(intersect, list(k1,k2))
# keep CpG sites with at least 4 reads mapped in each haplotype
BSobj_filt <- BSobj[keep_sites,]

dmlTest.sm.filt = DMLtest(BSobj_filt, group1=c("hap1.1","hap1.2","hap1.4","hap1.5","hap1.6"), 
                                  group2=c("hap2.1", "hap2.2","hap2.4","hap2.5","hap2.6"), smoothing=TRUE)

# call DMLs 
dmls = callDML(dmlTest.sm.filt, p.threshold=0.001) 
# call DMRs
# require at least 10 CpG sites to form a DMR + merge DMRs less than 70 bp apart
dmrs.minCG_10.filt = callDMR(dmlTest.sm.filt, p.threshold=0.001,minCG = 10,dis.merge = 70)
