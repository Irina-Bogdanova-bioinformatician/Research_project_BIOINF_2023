# Global gene expression of NSCLC TIL

Research project at the Bioinformatics Institute (2023).

### Objective
 Reproduction of the study of the transcription profile of tumor-infiltrating lymphocytes (TIL) described in the article ‘Transcriptional programs of neoantigen-specific TIL in anti-PD-1-treated lung cancers' (doi: 10.1038/s41586-021-03752-4).

### Tasks

**1. Data structure understanding and downloading**

Cell Ranger v3.1.0 was used by the authors to demultiplex the FASTQ reads, align them to the GRCh38 human transcriptome, and extract their cell and unique molecular identifier (UMI) barcodes. The output of this pipeline is a digital gene expression (DGE) matrix for each sample, which records the number of UMIs for each gene that are associated with each cell barcode. 
We downloaded preprocessed data from GEO (single-cell RNA-seq). Data structure of count matrix: one string is for one cell and one column is for one gene. The count matrix was read into an AnnData object, which holds many slots for annotations and different representations of the data. It also comes with its own HDF5-based file format: .h5ad. Scanpy toolkit (1.9.2) was used for analyzing single-cell gene expression. Scanpy is built jointly with anndata.

"MD01-019” - is a typical code of patient in this study. Any probe has alike prefix.
“MD01-005_tumor_5” -  is a typical code of probe, it has type of probe (tumor, normal, etc) and a number.

**2. Samples selection**

Six patients with and without major pathologic response were randomly selected:
MD01-019 (non-MPRs)
NY016-007 (non-MPRs)
NY016-014 (non-MPRs)
MD043-008 (MPR)
MD043-003 (MPR)
MD01-005 (MPR)

**3. QC** 

The quality of cells was then assessed based on (1) the number of genes detected per cell and (2) the proportion of mitochondrial gene/ribosomal gene counts. Low-quality cells were filtered if the number of detected genes was below 250 or above 3× the median absolute deviation away from the median gene number of all cells. Cells were filtered out if the proportion of mitochondrial gene counts was higher than 10% or the proportion of ribosomal genes was less than 10%. Mitochondrial genes (annotated with the prefix ‘MT-’), high abundance lincRNA genes, genes linked with poorly supported transcriptional models (annotated with the prefix ‘RP-’) and TCR (TR) genes (TRA/TRB/TRD/TRG, to avoid clonotype bias) were removed from further analysis. In addition, genes that were expressed in less than five cells were excluded.

The main problem we faced when working with big data was the lack of computing power.  
12 Gb of available RAM was not enough for processing the data from selected patient's sample. Extra memory in Google Colab was used after many crushes of notebook, however it still was not enough for full data selected. For this reason, the data were subsampled up to 1500 cells per sample before quality control and then combined into a single table for further analysis.

Before QC AnnData object we worked with was 61500 × 33538 (n_obs × n_vars) and after QC it was 48789 × 18309.

**4. Searching for variable genes, PCA, UMAP**

Scanpy was used to normalize the raw count data, identify highly variable features, scale features (normalization and algorithmization), and integrate samples. PCA (principal component analysis) was performed based on the most variable features identified. As there were more than ten thousand genes in analysis, we needed strong dimension reduction for later clusterization. PCA reveals the main axes of variation and denoises the data. Gene features associated with type I Interferon (IFN) response, immunoglobulin genes and specific mitochondrial related genes were excluded from clustering to avoid cell subsets driven by the above genes. Dimension reduction was done using the RunUMAP function. 

**5. Identification and annotation of T-cell clusters**

Leiden graph-clustering method was used for clusterization (community detection based on optimizing modularity) by Traag et al. (2018). 
Clusters were labeled based on the expression of the canonical immune cell markers.

**6. TIL expression profile assessment**

### Conclusions

Cells are succesfully clusterizated and annotated, expression profile is obtained.
