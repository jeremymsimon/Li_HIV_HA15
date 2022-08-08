# Load data into R/Seurat

library(fishpond)
library(Seurat)
library(tidyverse)
library(SingleCellExperiment)
library(Matrix)
library(glmGamPoi)
library(miQC)
library(SeuratWrappers)
library(flexmix)

# Set working directory
setwd("/proj/jmsimon/BARC/Jiang/CD4_HIVdelta6patched_scRNA_062822/Analysis_062922")


# Load S+A counts for RNA object
gex_q <- fishpond::loadFry(fryDir='/proj/jmsimon/BARC/Jiang/CD4_HIVdelta6patched_scRNA_062822/Jiang_CD4_HIVdelta6_GEX_quant_crlike',outputFormat="scRNA")
rna <- counts(gex_q)

# Load counts for HTO object
hto_q <- fishpond::loadFry(fryDir='/proj/jmsimon/BARC/Jiang/CD4_HIVdelta6patched_scRNA_062822/Jiang_CD4_HIVdelta6_HTO_quant_crlike')
hto <- counts(hto_q)

# Subset since we know we only have 9 HTOs, and it is the first 9 here
# If we don't do this, there will be false doublets called downstream!
hto <- hto[1:9,]

# As per various discussions, there is a lookup table internal to CellRanger that is needed to "translate" HTO feature barcodes such that they match the RNA feature barcodes
# https://www.biostars.org/p/9506747/
# This lookup table can be obtained from: 
# https://github.com/10XGenomics/cellranger/raw/master/lib/python/cellranger/barcodes/translation/3M-february-2018.txt.gz
# More here:
# https://kb.10xgenomics.com/hc/en-us/articles/360031133451-Why-is-there-a-discrepancy-in-the-3M-february-2018-txt-barcode-whitelist-



map <- read.table("/proj/jmsimon/BARC/Jiang/CD4_HIVdelta6_scRNA_041922/Analysis_042022/3M-february-2018.txt")
rownames(map) <- map$V1

colnames(rna) <- map[colnames(rna),2]
length(intersect(colnames(rna), colnames(hto)))
# [1] 9951





# Find cell barcodes in common between two assays
common.cells <- intersect(colnames(rna), colnames(hto))

# Now collapse transcripts to genes for the RNA data
tx2gene <- read.table("/proj/jmsimon/BARC/Jiang/CD4_HIVdelta6_scRNA_041922/Analysis_042022/gencode-v36_withExtraMito_HIV_p379delta6_tx2gene.tsv",header=F,sep="\t",col.names=c("tx","gene"))
exp.txId <- rownames(rna)
exp.geneId <- as.vector(tx2gene$gene[match(exp.txId, tx2gene$tx)])
exp.tx.grp <- t(sparse.model.matrix(~ 0 + exp.geneId))
exp.summarized <- exp.tx.grp %*% rna
rownames(exp.summarized) <- rownames(exp.summarized) %>% str_replace_all(".+.geneId","")

# Create seurat object with RNA data subset to just the cells for which we also have HTO information
jiang.seurat <- CreateSeuratObject(exp.summarized[, which(colnames(exp.summarized) %in% common.cells)])

# Add HTO info as a separate assay, subset to just the cells for which we also have expression information
jiang.seurat[["HTO"]] <- CreateAssayObject(counts = hto[, which(colnames(hto) %in% common.cells)])

# Normalize and scale HTO data, then demultiplex samples
jiang.seurat <- NormalizeData(jiang.seurat, assay = "HTO", normalization.method = "CLR")
jiang.seurat <- ScaleData(jiang.seurat, assay = "HTO", verbose = F)
jiang.seurat <- HTODemux(jiang.seurat, assay = "HTO", positive.quantile = 0.99)

# Plot signal for each HTO detected
Idents(jiang.seurat) <- "HTO_maxID"
pdf("Jiang_10X_HTO_alevinfry_HTO_ridgePlot.pdf")
RidgePlot(jiang.seurat, assay = "HTO", features = rownames(jiang.seurat[["HTO"]]), ncol = 3)
dev.off()

# Run PCA and UMAP to get an embedding representing just the HTO signal

VariableFeatures(jiang.seurat,assay="HTO") <- rownames(jiang.seurat[["HTO"]]@counts)
jiang.seurat <- RunPCA(jiang.seurat, reduction.name = "hto.pca", reduction.key = "HPC_", verbose = F, approx=FALSE, assay="HTO")
jiang.seurat <- RunUMAP(jiang.seurat, reduction = "hto.pca", dims = 1:9, reduction.name = "hto.umap", reduction.key = "HUMAP_", verbose = F, assay="HTO")

pdf("Jiang_10X_HTO_alevinfry_HTO_maxID_UMAP.pdf")
DimPlot(jiang.seurat, reduction = "hto.umap", label = T, group.by = "HTO_maxID")
dev.off()

pdf("Jiang_10X_HTO_alevinfry_HTO_hashID_UMAP.pdf")
DimPlot(jiang.seurat, reduction = "hto.umap", label = T, group.by = "hash.ID")
dev.off()

# Filter to just HTO singlets
jiang.seurat <- subset(jiang.seurat, subset = hash.ID != "Doublet")

# Filter to focus only on HA15 and negative control samples
# Added on 5/5/22
jiang.seurat <- subset(jiang.seurat, subset = hash.ID %in% c("Hashtag-1","Hashtag-2","Hashtag-3","Hashtag-4","Hashtag-5","Hashtag-6"))


# Run MiQC to remove high-mito cells
DefaultAssay(jiang.seurat) <- "RNA"
jiang.seurat <- PercentageFeatureSet(jiang.seurat, pattern = "^MT-", col.name = "percent.mt")
jiang.seurat <- RunMiQC(jiang.seurat, percent.mt = "percent.mt", nFeature_RNA = "nFeature_RNA", posterior.cutoff = 0.75, model.slot = "flexmix_model")
jiang.seurat <- subset(jiang.seurat, miQC.keep == "keep")

# Additionally filter low-count cells
jiang.seurat <- subset(jiang.seurat, subset = nCount_RNA > 2000 & nFeature_RNA > 1000)

dim(jiang.seurat)
[1] 59437  5312

# Add metadata for each sample
# Edited 5/5/22

annot = as.data.frame(cbind("HTO" = sort(unique(as.character(jiang.seurat$hash.ID))), "Sample" = c("NegControl_Rep1","NegControl_Rep2","NegControl_Rep3","HA15_Rep1","HA15_Rep2","HA15_Rep3")))
rownames(annot) = annot$HTO
annot$Treatment = cbind(str_replace_all(annot$Sample,"(.+)_.+","\\1"))
annot$Replicate = cbind(str_replace_all(annot$Sample,".+_(.+)","\\1"))

jiang.seurat <- AddMetaData(jiang.seurat,metadata=as.character(annot[as.character(jiang.seurat$hash.ID),2]),col.name="Sample")
jiang.seurat <- AddMetaData(jiang.seurat,metadata=as.character(annot[as.character(jiang.seurat$hash.ID),3]),col.name="Treatment")
jiang.seurat <- AddMetaData(jiang.seurat,metadata=as.character(annot[as.character(jiang.seurat$hash.ID),4]),col.name="Replicate")



# Split by HTO by biological sample
jiang.list = SplitObject(jiang.seurat, split.by="HTO_maxID")

# Get dims of each HTO
lapply(X = jiang.list, dim)

$`Hashtag-1`
[1] 59437   829

$`Hashtag-3`
[1] 59437   789

$`Hashtag-4`
[1] 59437  1316

$`Hashtag-2`
[1] 59437  1013

$`Hashtag-6`
[1] 59437   652

$`Hashtag-5`
[1] 59437   713



# Loop through list, run SCTransform as usual
jiang.list <- lapply(X = jiang.list, FUN = SCTransform, vst.flavor = "v2", return.only.var.genes = F, vars.to.regress = "percent.mt")
jiang.list <- lapply(X = jiang.list, FUN = FindVariableFeatures, selection.method = "vst", nfeatures = 3000, assay="SCT")

# Perform integration
options(future.globals.maxSize=5242880000)
all.features <- SelectIntegrationFeatures(object.list = jiang.list, nfeatures = 5000)
jiang.list <- PrepSCTIntegration(object.list = jiang.list, anchor.features = all.features)
all.anchors <- FindIntegrationAnchors(object.list = jiang.list, normalization.method = "SCT", anchor.features = all.features)
jiang.integrated <- IntegrateData(anchorset = all.anchors, normalization.method = "SCT")


#Save image at point of object integration
saveRDS(jiang.integrated,"Jiang_10X_HTO_alevinfry_SCT_062922_integrated.rds")


# Run PCA, then determine how many PCs are informative
jiang.integrated <- RunPCA(jiang.integrated, verbose = FALSE, npcs = 300)
ElbowPlot(jiang.integrated, ndims = 300, reduction = "pca")

# Run UMAP, and identify cell clusters
jiang.integrated <- RunUMAP(jiang.integrated, dims = 1:30)
jiang.integrated <- FindNeighbors(jiang.integrated, dims = 1:30, verbose = FALSE)
jiang.integrated <- FindClusters(jiang.integrated, verbose = FALSE, resolution = 0.25, algorithm=2)

# Plot UMAP labeled by clusters
pdf("Jiang_10X_HTO_alevinfry_SCT_062922_integrated_clustered_res0.25_umap.pdf")
DimPlot(jiang.integrated,reduction = "umap", label = TRUE)
DimPlot(jiang.integrated,reduction = "umap", group.by = "Sample")
DimPlot(jiang.integrated,reduction = "umap", group.by = "Treatment")
DimPlot(jiang.integrated,reduction = "umap", group.by = "Replicate")
dev.off()

saveRDS(jiang.integrated,"Jiang_10X_HTO_alevinfry_SCT_062922_integrated_clustered_res0.25.rds")


# Identify markers of each cluster
library(data.table)
library(presto)

vargenes<- presto::wilcoxauc(jiang.integrated, 'seurat_clusters', seurat_assay = 'SCT')
top_vargenes = top_markers(vargenes, n = 50, auc_min = 0.5, pct_in_min = 20, pct_out_max = 20)
write.table(top_vargenes,"Jiang_10X_HTO_alevinfry_SCT_062922_integrated_clustered_res0.25_prestoMarkers.txt",row.names=F,quote=F,sep="\t")

all_markers<- top_vargenes %>%
	select(-rank) %>% 
	unclass() %>% 
	stack() %>%
	pull(values) %>%
	unique() %>%
	.[!is.na(.)]



# The above list is too long for a dotplot, so filter further
top_vargenes = top_markers(vargenes, n = 10, auc_min = 0.5, pct_in_min = 20, pct_out_max = 20)
all_markers<- top_vargenes %>%
	select(-rank) %>% 
	unclass() %>% 
	stack() %>%
	pull(values) %>%
	unique() %>%
	.[!is.na(.)]

pdf("Jiang_10X_HTO_alevinfry_SCT_062922_integrated_clustered_res0.25_prestoMarkers_dotplot.pdf",width=30,height=6)
DotPlot(jiang.integrated,features=rev(all_markers),assay="SCT",cols = c("blue","red")) + RotatedAxis()
dev.off()


# Plot HIV genes
DefaultAssay(jiang.integrated) <- "SCT"

pdf("Jiang_10X_HTO_alevinfry_SCT_062922_integrated_clustered_res0.25_TatRev_FeaturePlot.pdf")
FeaturePlot(jiang.integrated,slot="data",features=c("Tat","Rev"),order=T)
dev.off()





# Tabulate proportions per cluster

proptable = as.data.frame(table(jiang.integrated$seurat_clusters, jiang.integrated$Sample)) %>% 
	as_tibble() %>%
	dplyr::rename("Sample" = Var2) %>%
	left_join(as_tibble(as.data.frame(table(jiang.integrated$Sample))) %>% dplyr::rename("Sample" = Var1), by="Sample") %>%
	mutate(Proportion = Freq.x / Freq.y) %>%
	select(-Freq.x,-Freq.y) %>%
	pivot_wider(names_from=Sample,values_from=Proportion) %>%
	dplyr::rename("Cluster" = Var1)


write.table(proptable,"Jiang_10X_HTO_alevinfry_SCT_062922_integrated_clustered_res0.25_proportions_replicates.txt",row.names=F,quote=FALSE,sep="\t")

# boxplot
pdf("Jiang_10X_HTO_alevinfry_SCT_062922_integrated_clustered_res0.25_proportions_replicates_boxplots.pdf",width=12,height=8)
proptable %>%
	pivot_longer(cols=!Cluster,names_to="Sample",values_to="Proportion") %>%
	tidyr::extract(Sample,c("Treatment","Replicate"),"(^.+)_(.+)") %>%
	mutate(Treatment = fct_relevel(Treatment,c("NegControl","HA15"))) %>%
    ggplot(aes(fill = Treatment, x = Treatment, y=Proportion)) +
    geom_boxplot(outlier.shape = NA) +
    geom_point(pch = 21, position = position_jitterdodge(jitter.width=0.2)) +
    facet_wrap(~Cluster,scales="free") +
    scale_fill_manual(values=c("#f5940c","#0fb7fa")) +
    xlab("")
dev.off()


# DA testing with DESeq2
library(DESeq2)

proptable = matrix(table(jiang.integrated$seurat_clusters, jiang.integrated$Sample),ncol=6)
rownames(proptable) = paste0("Cluster",seq(0,nrow(proptable)-1))
colnames(proptable) = colnames(table(jiang.integrated$seurat_clusters, jiang.integrated$Sample))


coldata = as.data.frame(cbind("Sample" = colnames(proptable), "Treatment" = str_replace_all(colnames(proptable),"(.+)_.+","\\1")))
rownames(coldata) = coldata$Sample
coldata$Sample <- factor(coldata$Sample)
coldata$Treatment <- factor(coldata$Treatment)

dds <- DESeqDataSetFromMatrix(countData = proptable, colData = coldata, design = ~ Treatment)
dds <- DESeq(dds,fitType="mean")
res <- results(dds, contrast=c("Treatment","HA15","NegControl"))
res <- res[order(res$padj), ]

rownames_to_column(as.data.frame(res),var="Cluster") %>% 
	as_tibble() %>% 
	filter(padj < 0.05)

  Cluster  baseMean log2FoldChange lfcSE  stat  pvalue   padj
  <chr>       <dbl>          <dbl> <dbl> <dbl>   <dbl>  <dbl>
1 Cluster5     22.0         -0.834 0.285 -2.92 0.00345 0.0276











# Now get Rev/Tat+ cells
# Use functions posted here:
# https://github.com/satijalab/seurat/issues/371


PrctCellExpringGene <- function(object, genes, group.by = "all"){
    if(group.by == "all"){
        prct = unlist(lapply(genes,calc_helper, object=object))
        result = data.frame(Markers = genes, Cell_proportion = prct)
        return(result)
    }

    else{        
        list = SplitObject(object, group.by)
        factors = names(list)

        results = lapply(list, PrctCellExpringGene, genes=genes)
        for(i in 1:length(factors)){
        results[[i]]$Feature = factors[i]
        }
        combined = do.call("rbind", results)
        return(combined)
    }
}

calc_helper <- function(object,genes){
    counts = object[['RNA']]@counts
    ncells = ncol(counts)
    if(genes %in% row.names(counts)){
    sum(counts[genes,]>0)/ncells
    }else{return(NA)}
}

# First compute and plot percent of cells expressing these two genes
table = c()
for (i in levels(jiang.integrated$seurat_clusters)) {
	clust.obj = subset(jiang.integrated,subset = seurat_clusters == i)
	tab = PrctCellExpringGene(clust.obj, genes =c("Tat","Rev"), group.by = "Sample")
	tab = tab %>%
		as_tibble() %>%
		add_column(Cluster = i)
	table = bind_rows(table,tab)
}

pdf("Jiang_10X_HTO_alevinfry_SCT_062922_integrated_clustered_res0.25_proportions_RevTat_expressing_boxplots.pdf",width=12,height=8)
table %>%
	separate(Feature,c("Treatment","Replicate"),sep="_",remove=F) %>%
	ggplot(aes(x=Treatment,y=Cell_proportion,fill=Treatment)) +
	geom_boxplot(outlier.shape=NA) + 
	geom_jitter(width=0.25) +
	facet_grid(Markers~Cluster)
dev.off()


# Now create a subsetable switch for expression of Rev/Tat

hivPos = colnames(jiang.integrated) %in% WhichCells(jiang.integrated, slot = 'counts', expression = Rev > 0 | Tat > 0 )
hivPos[hivPos==TRUE] <- "HIVpos"
hivPos[hivPos==FALSE] <- "HIVneg"

jiang.integrated = AddMetaData(jiang.integrated,metadata=hivPos,col.name="HIV")

# Check UMAP plot for HIV status label
DimPlot(jiang.integrated,reduction = "umap", label = TRUE, group.by = "HIV")

saveRDS(jiang.integrated,"Jiang_10X_HTO_alevinfry_SCT_062922_integrated_clustered_res0.25_withHIVlabels.rds")





# Compare HIV+ cell counts across the two treatment groups, per cluster


#jiang.integrated = readRDS("Jiang_10X_HTO_alevinfry_SCT_062922_integrated_clustered_res0.25_withHIVlabels.rds")

table = c()
for(i in levels(jiang.integrated$seurat_clusters)) {
	tab <- table(jiang.integrated$Sample[jiang.integrated$seurat_clusters==i],jiang.integrated$HIV[jiang.integrated$seurat_clusters==i]) %>% 
		as.data.frame() %>% 
		as_tibble() %>%
		dplyr::rename("Sample" = Var1, "HIV" = Var2) %>%
		add_column("Cluster" = paste0("Cluster",i))
	table = bind_rows(table,tab)
}

pdf("Jiang_10X_HTO_alevinfry_SCT_062922_integrated_clustered_res0.25_proportions_HIVposneg_byClusterTreatment_boxplots.pdf",width=12,height=8)
table %>%
	separate(Sample,c("Treatment","Replicate"),sep="_",remove=F) %>%
	group_by(Sample) %>%
	summarize(TotalN = sum(Freq)) %>%
	full_join(table %>% separate(Sample,c("Treatment","Replicate"),sep="_",remove=F),by=c("Sample")) %>%
	mutate(Ratio = Freq/TotalN) %>%
	ggplot(aes(x=Treatment,y=Ratio,fill=Treatment)) +
	geom_boxplot(outlier.shape=NA) + 
	geom_jitter(width=0.25) +
	facet_grid(HIV~Cluster,scales="free")
dev.off()


pdf("Jiang_10X_HTO_alevinfry_SCT_062922_integrated_clustered_res0.25_proportions_HIVposneg_byTreatment_boxplots.pdf")
table %>%
	separate(Sample,c("Treatment","Replicate"),sep="_",remove=F) %>%
	group_by(Sample) %>%
	summarize(TotalN = sum(Freq)) %>%
	full_join(table %>% separate(Sample,c("Treatment","Replicate"),sep="_",remove=F),by=c("Sample")) %>%
	mutate(Ratio = Freq/TotalN) %>%
	group_by(Sample,Treatment,HIV) %>%
	summarize(sum = sum(Ratio)) %>%
	ggplot(aes(x=Treatment,y=sum,fill=Treatment)) +
	geom_boxplot(outlier.shape=NA) + 
	geom_jitter(width=0.25) +
	facet_grid(~HIV,scales="free")
dev.off()











# Differential expression analysis on HA15 vs NegControl using Distinct, all cells regardless of HIV

library(distinct)
library(scuttle)
library(scater)
library(limma)
library(muscat)
library(gprofiler2)

# Convert to SingleCellExperiment, use raw counts
jiang.sce = as.SingleCellExperiment(jiang.integrated,assay="RNA")


# Compute log2 CPM consistent with distinct's code
jiang.sce = computeLibraryFactors(jiang.sce)
jiang.sce <- logNormCounts(jiang.sce)

# Compute plain CPM for log2FC calculations
cpm(jiang.sce) = calculateCPM(jiang.sce,assay.type="counts")

(jiang.sce <- prepSCE(jiang.sce,
	kid = "seurat_clusters",
	gid = "Treatment",
	sid = "Sample",
	drop=T))

colData(jiang.sce)

samples = jiang.sce@metadata$experiment_info$sample_id
group = jiang.sce@metadata$experiment_info$group_id
design = model.matrix(~group)
# rownames of the design must indicate sample ids:
rownames(design) = samples

set.seed(61217)
res = distinct_test(x = jiang.sce, 
                    name_assays_expression = "logcounts",
                    name_cluster = "cluster_id",
                    name_sample = "sample_id",
                    design = design,
                    column_to_test = 2,
                    min_non_zero_cells = 50,
                    n_cores = 2)
                    
res = log2_FC(res = res,
              x = jiang.sce, 
              name_assays_expression = "cpm",
              name_group = "group_id",
              name_cluster = "cluster_id")


prefix = "Jiang_10X_HTO_alevinfry_SCT_062922_integrated_clustered_res0.25_DEG_distinct_log2_cluster"
for(i in levels(jiang.sce$cluster_id)) {
	tr = top_results(res,global=F,cluster=as.character(i),significance=0.05,sort_by="p_adj.glb")
	if(length(tr[[1]])>0) {	
		file = paste0(prefix,i,".txt")
		print(file)
		write.table(tr,file,quote=F,sep="\t",row.names=F)
	}
}

save.image("Jiang_10X_HTO_alevinfry_SCT_062922_integrated_clustered_res0.25_DEG_distinct_log2.Rdata")




# Run gProfiler
for (i in as.numeric(levels(jiang.integrated$seurat_clusters))){		
	# Instantiate all variables to make sure no carry-overs from previous clusters
	names = 0
	diff = 0
	sig.up = 0
	sig.dn = 0
	gp = 0
	gp.flat = 0

	print(paste0("Working on cluster ",i))

	if(file.exists(paste0("Jiang_10X_HTO_alevinfry_SCT_062922_integrated_clustered_res0.25_DEG_distinct_log2_cluster",i,".txt"))) {
		diff = read.table(paste0("Jiang_10X_HTO_alevinfry_SCT_062922_integrated_clustered_res0.25_DEG_distinct_log2_cluster",i,".txt"),header=T,sep="\t",row.names=1)

		# Separate up and downregulated
		
		sig.up = rownames(diff[diff$p_adj.glb<0.05 & diff$log2FC_HA15.NegControl > 0 & !is.na(diff$log2FC_HA15.NegControl),])
		sig.dn = rownames(diff[diff$p_adj.glb<0.05 & diff$log2FC_HA15.NegControl < 0 & !is.na(diff$log2FC_HA15.NegControl),])

		if(length(sig.up)>=5 & length(sig.dn)>=5) {
			q = list(sig.up, sig.dn)
			names(q) = c(paste0("Cluster",i,"_UP"),paste0("Cluster",i,"_DN"))
			gp = gost(q,significant=TRUE,organism = "hsapiens",evcodes=TRUE,correction_method="fdr",sources=c("GO:BP","GO:MF","GO:CC","KEGG","REAC"),multi_query=FALSE)
			gp.flat = as_tibble(gp$result) %>%
				select(-parents,-evidence_codes)	
			write.table(gp.flat,file=paste0("Jiang_10X_HTO_alevinfry_SCT_062922_integrated_clustered_res0.25_DEG_distinct_log2_cluster",i,"_gProfiler_051022.txt"),quote=F,sep="\t",col.names=NA)
		} else if(length(sig.up)>=5) {
			q = list(sig.up)
			names(q) = paste0("Cluster",i,"_UP")
			gp = gost(q,significant=TRUE,organism = "hsapiens",evcodes=TRUE,correction_method="fdr",sources=c("GO:BP","GO:MF","GO:CC","KEGG","REAC"),multi_query=FALSE)
			gp.flat = as_tibble(gp$result) %>%
				select(-parents,-evidence_codes)	
			write.table(gp.flat,file=paste0("Jiang_10X_HTO_alevinfry_SCT_062922_integrated_clustered_res0.25_DEG_distinct_log2_cluster",i,"_gProfiler_051022.txt"),quote=F,sep="\t",col.names=NA)
		} else if(length(sig.dn)>=5) {
			q = list(sig.dn)
			names(q) = paste0("Cluster",i,"_DN")
			gp = gost(q,significant=TRUE,organism = "hsapiens",evcodes=TRUE,correction_method="fdr",sources=c("GO:BP","GO:MF","GO:CC","KEGG","REAC"),multi_query=FALSE)
			gp.flat = as_tibble(gp$result) %>%
				select(-parents,-evidence_codes)	
			write.table(gp.flat,file=paste0("Jiang_10X_HTO_alevinfry_SCT_062922_integrated_clustered_res0.25_DEG_distinct_log2_cluster",i,"_gProfiler_051022.txt"),quote=F,sep="\t",col.names=NA)
		} else {
			# skip
		}
	}
}








































# Differential expression analysis on HA15 vs NegControl using Distinct, HIV+ cells only

jiang.hiv <- subset(jiang.integrated,subset=HIV=="HIVpos")

# Convert to SingleCellExperiment, use raw counts
jiang.sce = as.SingleCellExperiment(jiang.hiv,assay="RNA")


# Compute log2 CPM consistent with distinct's code
jiang.sce = computeLibraryFactors(jiang.sce)
jiang.sce <- logNormCounts(jiang.sce)

# Compute plain CPM for log2FC calculations
cpm(jiang.sce) = calculateCPM(jiang.sce,assay.type="counts")

(jiang.sce <- prepSCE(jiang.sce,
	kid = "seurat_clusters",
	gid = "Treatment",
	sid = "Sample",
	drop=T))

colData(jiang.sce)

samples = jiang.sce@metadata$experiment_info$sample_id
group = jiang.sce@metadata$experiment_info$group_id
design = model.matrix(~group)
# rownames of the design must indicate sample ids:
rownames(design) = samples

set.seed(61217)
res = distinct_test(x = jiang.sce, 
                    name_assays_expression = "logcounts",
                    name_cluster = "cluster_id",
                    name_sample = "sample_id",
                    design = design,
                    column_to_test = 2,
                    min_non_zero_cells = 50,
                    n_cores = 2)
                    
res = log2_FC(res = res,
              x = jiang.sce, 
              name_assays_expression = "cpm",
              name_group = "group_id",
              name_cluster = "cluster_id")


prefix = "Jiang_10X_HTO_alevinfry_SCT_062922_integrated_clustered_res0.25_HIVpos_DEG_distinct_log2_cluster"
for(i in levels(jiang.sce$cluster_id)) {
	tr = top_results(res,global=F,cluster=as.character(i),significance=0.05,sort_by="p_adj.glb")
	if(length(tr[[1]])>0) {	
		file = paste0(prefix,i,".txt")
		print(file)
		write.table(tr,file,quote=F,sep="\t",row.names=F)
	}
}

save.image("Jiang_10X_HTO_alevinfry_SCT_062922_integrated_clustered_res0.25_HIVpos_DEG_distinct_log2.Rdata")




# Run gProfiler
for (i in as.numeric(levels(jiang.integrated$seurat_clusters))){		
	# Instantiate all variables to make sure no carry-overs from previous clusters
	names = 0
	diff = 0
	sig.up = 0
	sig.dn = 0
	gp = 0
	gp.flat = 0

	print(paste0("Working on cluster ",i))

	if(file.exists(paste0("Jiang_10X_HTO_alevinfry_SCT_062922_integrated_clustered_res0.25_HIVpos_DEG_distinct_log2_cluster",i,".txt"))) {
		diff = read.table(paste0("Jiang_10X_HTO_alevinfry_SCT_062922_integrated_clustered_res0.25_HIVpos_DEG_distinct_log2_cluster",i,".txt"),header=T,sep="\t",row.names=1)

		# Separate up and downregulated
		
		sig.up = rownames(diff[diff$p_adj.glb<0.05 & diff$log2FC_HA15.NegControl > 0 & !is.na(diff$log2FC_HA15.NegControl),])
		sig.dn = rownames(diff[diff$p_adj.glb<0.05 & diff$log2FC_HA15.NegControl < 0 & !is.na(diff$log2FC_HA15.NegControl),])

		if(length(sig.up)>=5 & length(sig.dn)>=5) {
			q = list(sig.up, sig.dn)
			names(q) = c(paste0("Cluster",i,"_UP"),paste0("Cluster",i,"_DN"))
			gp = gost(q,significant=TRUE,organism = "hsapiens",evcodes=TRUE,correction_method="fdr",sources=c("GO:BP","GO:MF","GO:CC","KEGG","REAC"),multi_query=FALSE)
			gp.flat = as_tibble(gp$result) %>%
				select(-parents,-evidence_codes)	
			write.table(gp.flat,file=paste0("Jiang_10X_HTO_alevinfry_SCT_062922_integrated_clustered_res0.25_HIVpos_DEG_distinct_log2_cluster",i,"_gProfiler_051022.txt"),quote=F,sep="\t",col.names=NA)
		} else if(length(sig.up)>=5) {
			q = list(sig.up)
			names(q) = paste0("Cluster",i,"_UP")
			gp = gost(q,significant=TRUE,organism = "hsapiens",evcodes=TRUE,correction_method="fdr",sources=c("GO:BP","GO:MF","GO:CC","KEGG","REAC"),multi_query=FALSE)
			gp.flat = as_tibble(gp$result) %>%
				select(-parents,-evidence_codes)	
			write.table(gp.flat,file=paste0("Jiang_10X_HTO_alevinfry_SCT_062922_integrated_clustered_res0.25_HIVpos_DEG_distinct_log2_cluster",i,"_gProfiler_051022.txt"),quote=F,sep="\t",col.names=NA)
		} else if(length(sig.dn)>=5) {
			q = list(sig.dn)
			names(q) = paste0("Cluster",i,"_DN")
			gp = gost(q,significant=TRUE,organism = "hsapiens",evcodes=TRUE,correction_method="fdr",sources=c("GO:BP","GO:MF","GO:CC","KEGG","REAC"),multi_query=FALSE)
			gp.flat = as_tibble(gp$result) %>%
				select(-parents,-evidence_codes)	
			write.table(gp.flat,file=paste0("Jiang_10X_HTO_alevinfry_SCT_062922_integrated_clustered_res0.25_HIVpos_DEG_distinct_log2_cluster",i,"_gProfiler_051022.txt"),quote=F,sep="\t",col.names=NA)
		} else {
			# skip
		}
	}
}







































# Differential expression analysis on HA15 vs NegControl using Distinct, HIV- cells only

jiang.hiv <- subset(jiang.integrated,subset=HIV=="HIVneg")

# Convert to SingleCellExperiment, use raw counts
jiang.sce = as.SingleCellExperiment(jiang.hiv,assay="RNA")


# Compute log2 CPM consistent with distinct's code
jiang.sce = computeLibraryFactors(jiang.sce)
jiang.sce <- logNormCounts(jiang.sce)

# Compute plain CPM for log2FC calculations
cpm(jiang.sce) = calculateCPM(jiang.sce,assay.type="counts")

(jiang.sce <- prepSCE(jiang.sce,
	kid = "seurat_clusters",
	gid = "Treatment",
	sid = "Sample",
	drop=T))

colData(jiang.sce)

samples = jiang.sce@metadata$experiment_info$sample_id
group = jiang.sce@metadata$experiment_info$group_id
design = model.matrix(~group)
# rownames of the design must indicate sample ids:
rownames(design) = samples

set.seed(61217)
res = distinct_test(x = jiang.sce, 
                    name_assays_expression = "logcounts",
                    name_cluster = "cluster_id",
                    name_sample = "sample_id",
                    design = design,
                    column_to_test = 2,
                    min_non_zero_cells = 50,
                    n_cores = 2)
                    
res = log2_FC(res = res,
              x = jiang.sce, 
              name_assays_expression = "cpm",
              name_group = "group_id",
              name_cluster = "cluster_id")


prefix = "Jiang_10X_HTO_alevinfry_SCT_062922_integrated_clustered_res0.25_HIVneg_DEG_distinct_log2_cluster"
for(i in levels(jiang.sce$cluster_id)) {
	tr = top_results(res,global=F,cluster=as.character(i),significance=0.05,sort_by="p_adj.glb")
	if(length(tr[[1]])>0) {	
		file = paste0(prefix,i,".txt")
		print(file)
		write.table(tr,file,quote=F,sep="\t",row.names=F)
	}
}

save.image("Jiang_10X_HTO_alevinfry_SCT_062922_integrated_clustered_res0.25_HIVneg_DEG_distinct_log2.Rdata")




# Run gProfiler
for (i in as.numeric(levels(jiang.integrated$seurat_clusters))){		
	# Instantiate all variables to make sure no carry-overs from previous clusters
	names = 0
	diff = 0
	sig.up = 0
	sig.dn = 0
	gp = 0
	gp.flat = 0

	print(paste0("Working on cluster ",i))

	if(file.exists(paste0("Jiang_10X_HTO_alevinfry_SCT_062922_integrated_clustered_res0.25_HIVneg_DEG_distinct_log2_cluster",i,".txt"))) {
		diff = read.table(paste0("Jiang_10X_HTO_alevinfry_SCT_062922_integrated_clustered_res0.25_HIVneg_DEG_distinct_log2_cluster",i,".txt"),header=T,sep="\t",row.names=1)

		# Separate up and downregulated
		
		sig.up = rownames(diff[diff$p_adj.glb<0.05 & diff$log2FC_HA15.NegControl > 0 & !is.na(diff$log2FC_HA15.NegControl),])
		sig.dn = rownames(diff[diff$p_adj.glb<0.05 & diff$log2FC_HA15.NegControl < 0 & !is.na(diff$log2FC_HA15.NegControl),])

		if(length(sig.up)>=5 & length(sig.dn)>=5) {
			q = list(sig.up, sig.dn)
			names(q) = c(paste0("Cluster",i,"_UP"),paste0("Cluster",i,"_DN"))
			gp = gost(q,significant=TRUE,organism = "hsapiens",evcodes=TRUE,correction_method="fdr",sources=c("GO:BP","GO:MF","GO:CC","KEGG","REAC"),multi_query=FALSE)
			gp.flat = as_tibble(gp$result) %>%
				select(-parents,-evidence_codes)	
			write.table(gp.flat,file=paste0("Jiang_10X_HTO_alevinfry_SCT_062922_integrated_clustered_res0.25_HIVneg_DEG_distinct_log2_cluster",i,"_gProfiler_051022.txt"),quote=F,sep="\t",col.names=NA)
		} else if(length(sig.up)>=5) {
			q = list(sig.up)
			names(q) = paste0("Cluster",i,"_UP")
			gp = gost(q,significant=TRUE,organism = "hsapiens",evcodes=TRUE,correction_method="fdr",sources=c("GO:BP","GO:MF","GO:CC","KEGG","REAC"),multi_query=FALSE)
			gp.flat = as_tibble(gp$result) %>%
				select(-parents,-evidence_codes)	
			write.table(gp.flat,file=paste0("Jiang_10X_HTO_alevinfry_SCT_062922_integrated_clustered_res0.25_HIVneg_DEG_distinct_log2_cluster",i,"_gProfiler_051022.txt"),quote=F,sep="\t",col.names=NA)
		} else if(length(sig.dn)>=5) {
			q = list(sig.dn)
			names(q) = paste0("Cluster",i,"_DN")
			gp = gost(q,significant=TRUE,organism = "hsapiens",evcodes=TRUE,correction_method="fdr",sources=c("GO:BP","GO:MF","GO:CC","KEGG","REAC"),multi_query=FALSE)
			gp.flat = as_tibble(gp$result) %>%
				select(-parents,-evidence_codes)	
			write.table(gp.flat,file=paste0("Jiang_10X_HTO_alevinfry_SCT_062922_integrated_clustered_res0.25_HIVneg_DEG_distinct_log2_cluster",i,"_gProfiler_051022.txt"),quote=F,sep="\t",col.names=NA)
		} else {
			# skip
		}
	}
}

























# Explore expression of genes in Ferroptosis pathway, separated by Cluster, Treatment, and HIV status
library(circlize)
library(tidyHeatmap)
library(Seurat)
library(tidyverse)
library(distinct)
library(scuttle)
library(scater)
library(limma)
library(muscat)

#jiang.integrated = readRDS("Jiang_10X_HTO_alevinfry_SCT_062922_integrated_clustered_res0.25_withHIVlabels.rds")

ferro <- c("GCLC","GCLM","CP","ATG5","ACSL4","TFRC","VDAC3","FTL","TF","HMOX1","GSS","MAP1LC3A","SLC39A14","SLC11A2","LPCAT3","STEAP3","ACSL3","SAT1","SLC40A1","SLC39A8","MAP1LC3B","SAT2","TP53","SLC7A11","ACSL1","ALOX15","ACSL6","CYBB","VDAC2","GPX4","FTH1","SLC3A2","PCBP1","PRNP","FTMT","PCBP2","ACSL5","ATG7","MAP1LC3C","MAP1LC3B2","NCOA4")

# Convert to SingleCellExperiment, use raw counts
jiang.sce = as.SingleCellExperiment(jiang.integrated,assay="RNA")


# Compute log2 CPM consistent with distinct's code
jiang.sce = computeLibraryFactors(jiang.sce)
jiang.sce <- logNormCounts(jiang.sce)

# Now convert back to Seurat object format
jiang.cpm.seurat <- as.Seurat(jiang.sce, counts = "counts", data = "logcounts")




# Tidy data and plot as a heatmap
# This is slow but keeping it sequential for readability

tidy = rownames_to_column(as.data.frame(jiang.cpm.seurat@assays$RNA@data),var="Gene") %>% 
	as_tibble() %>%
	pivot_longer(cols=!Gene,names_to="Cell",values_to="CPM") %>%
	full_join(
		bind_rows(jiang.cpm.seurat$HIV) %>% 
			pivot_longer(cols=everything()) %>%
			dplyr::rename("Cell" = name, "HIV" = value), by="Cell"	
	) %>%
	full_join(
		bind_rows(jiang.cpm.seurat$seurat_clusters) %>% 
			pivot_longer(cols=everything()) %>%
			dplyr::rename("Cell" = name, "Cluster" = value), by="Cell"	
	) %>%
	full_join(
		bind_rows(jiang.cpm.seurat$Treatment) %>% 
			pivot_longer(cols=everything()) %>%
			dplyr::rename("Cell" = name, "Treatment" = value), by="Cell"	
	) %>%
	full_join(
		bind_rows(jiang.cpm.seurat$Sample) %>% 
			pivot_longer(cols=everything()) %>%
			dplyr::rename("Cell" = name, "Sample" = value), by="Cell"	
	) %>%
	filter(Gene %in% ferro)
	

# Set column ordering, just use first gene
col_order <- tidy %>%
	group_by(HIV) %>%
	arrange(Cluster,desc(Treatment)) %>%
	filter(Gene=="ACSL1") %>%
	pull(Cell)

pdf("Jiang_10X_HTO_alevinfry_SCT_062922_integrated_clustered_res0.25_ferroptosis_CPM_tidyHeatmap.pdf",width=12,height=8)
tidy %>%
	filter(!Gene %in% c("STEAP3","MAP1LC3C","FTMT","ALOX15","CYBB")) %>%		# Remove these genes as they are not expressed
	group_by(HIV) %>%
	tidyHeatmap::heatmap(Gene,Cell,CPM,.scale="row",clustering_distance_rows = function(x) as.dist(1-cor(t(x))),clustering_method_rows="complete",palette_value=colorRamp2(c(-1, 0, 1), c("white", "white", "red")),cluster_columns=F,column_order=col_order,use_raster=T,raster_quality=2) %>%
	add_tile(Cluster) %>%
	add_tile(Treatment)
dev.off()





# Export cluster identities and SCT-normalized data matrix for GEO

write.table(as.data.frame(jiang.integrated$seurat_clusters),"ClusterIdents.txt",quote=F,sep="\t",col.names=NA)
write.table(as.data.frame(jiang.integrated$Sample),"SampleIdents.txt",quote=F,sep="\t",col.names=NA)
write.table(jiang.integrated@assays$SCT@data,"SCTnormalized.txt",quote=F,sep="\t",col.names=NA)
