scRNA-seq data preprocessing for *seismic* input
================
Qiliang Lai
9/16/2024

## Introduction

The *seismic* framework takes the SingleCellExperiment object as input
while requiring a column of cell metadata to specify the analysis
granularity. Before the *seismic* analysis, the data should be
preprocessed to ensure the data quality and compatibility with the
*seismic* framework.

## Step1: Load the scRNA-seq data and metadata

Assume now we already have the scRNA-seq data and metadata in the form
of SingleCellExperiment object. We will follow some tips from
[single-cell best
practice](https://www.sc-best-practices.org/preprocessing_visualization/normalization.html).

``` r
#load packages
library(SingleCellExperiment)
library(dplyr)
library(magrittr)
library(seismicGWAS)
```

``` r
#load data
example_sce <- readRDS("data/expr/sample_expr/sample_expr.rds") #your data path

example_sce
```

    ## class: SingleCellExperiment 
    ## dim: 23341 10000 
    ## metadata(0):
    ## assays(2): counts logcounts
    ## rownames(23341): 0610005C13Rik 0610007C21Rik ... l7Rn6
    ##   zsGreen-transgene
    ## rowData names(0):
    ## colnames: NULL
    ## colData names(39): nReads orig.ident ... res.0.2 ident
    ## reducedDimNames(0):
    ## mainExpName: RNA
    ## altExpNames(0):

## Step2: Check the data quality

Then we can check the data quality by looking at the distribution of the
gene expression and the number of genes expressed in each cell. Some
widely used metric to evaluate the data quality includes: sequencing
depth, number of genes detected, mitochondrial gene percentage, etc.

``` r
#total counts
example_sce$tot_counts <- colSums(assay(example_sce, "counts"))

#total mitrochondrial counts
example_sce$mito_ratio <- colSums(assay(example_sce, "counts")[grep("^Mt-", rownames(example_sce)), ]) / example_sce$tot_counts

#number of genes detected
example_sce$detected_genes <- colSums(assay(example_sce, "counts") > 0)
```

Typically these mitochondrial genes will be named with “Mt-” prefix in
mouse and “MT-” prefix in human. Here in these data these genes are not
present. Thus we only filter cells based on the total counts and the
number of detected genes. In addition, in this dataset we have cell
ontology information, indicating cell type annotation, so we should also
remove cells without such information.

``` r
#filter by total counts and detected genes
example_sce <- example_sce[, example_sce$tot_counts > 2000 & example_sce$detected_genes > 2000]

#filter by cell ontology information
example_sce <- example_sce[,!is.na(example_sce$cell_ontology_id)]

example_sce
```

    ## class: SingleCellExperiment 
    ## dim: 23341 7346 
    ## metadata(0):
    ## assays(2): counts logcounts
    ## rownames(23341): 0610005C13Rik 0610007C21Rik ... l7Rn6
    ##   zsGreen-transgene
    ## rowData names(0):
    ## colnames: NULL
    ## colData names(42): nReads orig.ident ... mito_ratio detected_genes
    ## reducedDimNames(0):
    ## mainExpName: RNA
    ## altExpNames(0):

### Optional to do: Remove genes with low expression to alleviate memory and computation burden

Also we can remove genes with low expression, because these genes are
less likely to be informative for downstream analysis. Here the gene
that is expressed in fewer than 5 cells or have fewer than 10 total RNA
counts detected are removed.

``` r
#number of cells expressing each gene
rowData(example_sce)$num_cells <- rowSums(assay(example_sce ,"counts")>0)

#number of total counts deteced for each gene
rowData(example_sce)$num_counts <- rowSums(assay(example_sce ,"counts"))

example_sce <- example_sce[rowData(example_sce)$num_cells>=5 & rowData(example_sce)$num_counts>=10,]
```

## Step3: Normalize the data

As desribed in [single-cell best
practice](https://www.sc-best-practices.org/preprocessing_visualization/normalization.html),
[scran](https://bioconductor.org/packages/release/bioc/html/scran.html)
was recommended for normalization, which estimates the per-cell
normalization factor based on cell pooling.

``` r
#imprt package
library(scran)
```

    ## Loading required package: scuttle

``` r
#cell pooling (cluster)
cell_pooling <- quickCluster(example_sce,  assay.type = "counts") 

#factor calculation
size_factor <- calculateSumFactors(example_sce, cluster=cell_pooling, min.mean=0.1, assay.type="counts")

#normalize
example_sce <- logNormCounts(example_sce, size.factors = size_factor )
```

## Step4: Choose analysis granularity

The *seismic* framework requires a column of cell metadata to specify
the analysis granularity. In current data we would like to care not only
about cell types but also tissue-specific effects, so we combine them
together.

``` r
example_sce$cell_type <- ifelse(!is.na(example_sce$free_annotation), 
        paste0(example_sce$tissue,".",example_sce$free_annotation), paste0(example_sce$tissue,".",example_sce$cell_ontology_class))
```

## Step5: Run the *seismic* analysis

Then we go through the normal *seismic* analysis pipeline. To refer to
the tutorial, please visit
[here](https://github.com/ylaboratory/seismic/blob/gh_page/vignettes/seismicGWAS.md).
The same processed gene-level MAGMA z-score for type 2 diabetes are
used.

``` r
#calculate specificity score
sscore  <- calc_specificity(sce = example_sce, ct_label_col = "cell_type")
#map to human genes
sscore_hsa <- translate_gene_ids(sscore, from = "mmu_symbol")
#association
ct_association <- get_ct_trait_associations(sscore = sscore_hsa, magma = t2d_magma)
#head
head(ct_association, n = 3)
```

    ##                     cell_type       pvalue         FDR
    ##                        <char>        <num>       <num>
    ## 1:         Pancreas.beta cell 6.330263e-06 0.000481100
    ## 2: Pancreas.pancreatic A cell 6.691778e-05 0.002542876
    ## 3: Pancreas.pancreatic D cell 4.749954e-04 0.012033218
