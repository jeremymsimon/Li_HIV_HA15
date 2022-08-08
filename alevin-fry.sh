# Largely following the tutorial here:
# https://combine-lab.github.io/alevin-fry-tutorials/2021/af-feature-bc/
# Note: example data uses chromium v2 chemistry

# Note 2: 
# Some of this was already done in the earlier iteration using the outdated HIV reference

# Build splici index for quantification
# Include extra mito sequences as these are absent in the gencode reference
# Downloaded from: https://zenodo.org/record/5799568 on 2/24/22
# In R:

library(roe)

setwd("/proj/jmsimon/BARC/Jiang/hg38-HIVp379d6")

make_splici_txome(
  gtf_path = "/proj/jmsimon/genomeAnnotation/gencode-v36_HIV_p379delta6_transcripts.gtf",
  genome_path = "/proj/jmsimon/genomeAnnotation/hg38-HIVp379d6.fasta",
  read_length = 90,
  flank_trim_length = 5,
  output_dir = "/proj/jmsimon/BARC/Jiang/hg38-HIVp379d6/",
  file_name_prefix = "gencode-v36_HIV_p379delta6_transcriptome_splici",
  extra_spliced = "/proj/jmsimon/genomeAnnotation/homo_sapiens_mito_seqs_for_splici.fasta",
  dedup_seqs = FALSE
)

# Then build salmon index as usual
ssub -n 12 --mem 100g --wrap=\"salmon index -t /proj/jmsimon/BARC/Jiang/hg38-HIVp379d6/gencode-v36_HIV_p379delta6_transcriptome_splici_fl85.fa -p 12 -i /proj/jmsimon/BARC/Jiang/hg38-HIVp379d6/gencode-v36_HIV_p379delta6_transcriptome_splici_fl85_idx\"


# Gather and index antibody sequences
wget --content-disposition  -nv https://ftp.ncbi.nlm.nih.gov/geo/series/GSE128nnn/GSE128639/suppl/GSE128639_MNC_HTO_Barcodes.csv.gz &&

gunzip -c GSE128639_MNC_HTO_Barcodes.csv.gz | awk -F "," '{print $1"\t"$4}' | sed 's/Hashtag /Hashtag_/g' | tail -n +2 > hto.tsv

salmon index -t hto.tsv -i hto_index --features -k7


# Now map the reads with alevin, RNA first
ssub -n 4 --mem 50g --wrap=\"salmon alevin -l ISR -i /proj/jmsimon/BARC/Jiang/hg38-HIVp379d6/gencode-v36_HIV_p379delta6_transcriptome_splici_fl85_idx \
  -1 /proj/gjianglab/HTSF/220404_UNC41-A00434_0452_AH323FDRX2/H323FDRX2/*GEX*R1*.fastq.gz \
  -2 /proj/gjianglab/HTSF/220404_UNC41-A00434_0452_AH323FDRX2/H323FDRX2/*GEX*R2*.fastq.gz \
  --chromiumV3 \
  -o /proj/jmsimon/BARC/Jiang/CD4_HIVdelta6_scRNA_062822/Jiang_CD4_HIVdelta6_GEX \
  -p 4 \
  --sketch\"
  
  
# Now HTO
# Some additional info about geometries here:
# https://github.com/ATpoint/sc_preprocess
ssub --mem 10g --wrap=\"salmon alevin -l ISR -i /proj/jmsimon/BARC/Jiang/hg38-HIVp379d6/hto_index \
  -1 /proj/gjianglab/HTSF/220404_UNC41-A00434_0452_AH323FDRX2/H323FDRX2/*CSP*R1*.fastq.gz \
  -2 /proj/gjianglab/HTSF/220404_UNC41-A00434_0452_AH323FDRX2/H323FDRX2/*CSP*R2*.fastq.gz \
  --read-geometry 2[11-25] --bc-geometry 1[1-16] --umi-geometry 1[17-28] \
  -o /proj/jmsimon/BARC/Jiang/CD4_HIVdelta6_scRNA_041922/Jiang_CD4_HIVdelta6_HTO --sketch\"



# Now quantify with alevin-fry

# Quantify RNA
ssub --mem 50g --wrap=\"alevin-fry generate-permit-list -d fw -i /proj/jmsimon/BARC/Jiang/CD4_HIVdelta6_scRNA_062822/Jiang_CD4_HIVdelta6_GEX -o /proj/jmsimon/BARC/Jiang/CD4_HIVdelta6_scRNA_062822/Jiang_CD4_HIVdelta6_GEX_quant -k\"
ssub --mem 50g --wrap=\"alevin-fry collate -r /proj/jmsimon/BARC/Jiang/CD4_HIVdelta6_scRNA_062822/Jiang_CD4_HIVdelta6_GEX -i /proj/jmsimon/BARC/Jiang/CD4_HIVdelta6_scRNA_062822/Jiang_CD4_HIVdelta6_GEX_quant -t 16\"
ssub --mem 50g --wrap=\"alevin-fry quant -m /proj/jmsimon/BARC/Jiang/hg38-HIVp379d6/gencode-v36_HIV_p379delta6_transcriptome_splici_fl85_t2g_3col.tsv -i /proj/jmsimon/BARC/Jiang/CD4_HIVdelta6_scRNA_062822/Jiang_CD4_HIVdelta6_GEX_quant -o /proj/jmsimon/BARC/Jiang/CD4_HIVdelta6_scRNA_062822/Jiang_CD4_HIVdelta6_GEX_quant_crlike -r cr-like -t 16 --use-mtx\"

# Now ADT/HTO
# But first we need fake tx2gene tables

awk '{print $1"\t"$1;}' hto.tsv > t2g_hto.tsv

# Quantify HTO
ssub --mem 50g --wrap=\"alevin-fry generate-permit-list -d fw -i /proj/jmsimon/BARC/Jiang/CD4_HIVdelta6_scRNA_041922/Jiang_CD4_HIVdelta6_HTO -o /proj/jmsimon/BARC/Jiang/CD4_HIVdelta6_scRNA_041922/Jiang_CD4_HIVdelta6_HTO_quant -k\"
ssub --mem 50g --wrap=\"alevin-fry collate -r /proj/jmsimon/BARC/Jiang/CD4_HIVdelta6_scRNA_041922/Jiang_CD4_HIVdelta6_HTO -i /proj/jmsimon/BARC/Jiang/CD4_HIVdelta6_scRNA_041922/Jiang_CD4_HIVdelta6_HTO_quant\"
ssub --mem 50g --wrap=\"alevin-fry quant -m /proj/jmsimon/BARC/Jiang/hg38-HIVp379d6/t2g_hto.tsv -i /proj/jmsimon/BARC/Jiang/CD4_HIVdelta6_scRNA_041922/Jiang_CD4_HIVdelta6_HTO_quant -o /proj/jmsimon/BARC/Jiang/CD4_HIVdelta6_scRNA_041922/Jiang_CD4_HIVdelta6_HTO_quant_crlike -r cr-like --use-mtx\"
