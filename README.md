# Li _et al._ 2022 (_manuscript under review_)
This repository contains code to support the scRNA-seq data analysis portion of the manuscript Li _et al._ 2022 (_manuscript under review_)
## What you'll find here:
* `alevin-fry.sh`
  - This file contains `bash` and `R` code used to construct a `splici` reference genome, download HTO barcode sequences, map GEX and HTO barcodes with `alevin`, then quantify GEX and HTO libraries with `alevin-fry`
  
* `seurat.R`
  - This file contains `R` code used to import scRNA-seq data and analyze with `Seurat`

## What you'll find elsewhere:
* Gene expression omnibus (GEO) accession GSEXXXXXXX contains all raw FASTQ files, `alevin-fry` output files, the SCT-normalized data matrix, as well as cluster and sample identities for each filtered barcode. This accession also contains the `gtf` and `fasta` file for the custom HIV reference genome, which is treated here as a pseudo-chromosome concatenated with hg38/GENCODE v36
