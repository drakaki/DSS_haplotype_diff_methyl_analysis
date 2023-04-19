# DSS_haplotype_diff_methyl_analysis
Differential methylation analysis between haplotypes using Bioconductor's DSS package for dispersion shrinkage

Citation

Wu H, Wang C, Wu Z (2013). “A new shrinkage estimator for dispersion improves differential expression detection in RNA-seq data.” 
Biostatistics. doi: 10.1093/biostatistics/kxs033.

Required file format

bed files from BS-seq or long read sequencing technologies (e.g. PacBio) with DNA polymerase kinetics, allowing for 5mC detection.

bed files uploaded to the Rstudio environment should contain four columns:
- chromosome
- 5mC site coordinate
- number of reads mapped to site
- number of mapped reads that are methylated,
as explained in the DSS vignette
