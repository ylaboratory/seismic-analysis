GWAS preprocessing for *seismic* input
================
Qiliang Lai
9/13/2024

## Introduction

The *seismic* framework takes the output from
[MAGMA](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1004219)
as input. This tutorial will guide you through the preprocessing steps
necessary to convert raw GWAS summary statistics into the proper format
for *seismic* analysis. Before beginning, ensure you have installed the
MAGMA software and obtained the required files: - Gene location file -
SNP synonyms file (dbSNP) - Population-specific reference data These
files can be all downloaded from the [MAGMA
website](https://cncr.nl/research/magma/).

## Step 1: Clean up the GWAS summary statistics file

The minimum required information in your GWAS summary statistics file
includes: - Chromosome ID column - Base-pair position column - P value
column Additionally, we strongly recommend including SNP IDs (rsIDs) to
facilitate MAGMA analysis. Here’s an example of how your data might
look:

``` r
#load packages
library(magrittr)
library(dplyr)
library(rtracklayer) #for lift over
```

``` r
#load data
example_gwas <- read.table("data/gwas/sample_gwas/sample_gwas.txt", header = T) %>% as_tibble() #your data path

head(example_gwas)
```

    ## # A tibble: 6 × 7
    ##   Allele1 Allele2 Direction P.value  group Chromosome     BP
    ##   <chr>   <chr>   <chr>       <dbl>  <int>      <int>  <int>
    ## 1 t       c       ---        0.0541 184129          1 796338
    ## 2 a       g       +++        0.0315  46089          1 817341
    ## 3 t       c       +++        0.0159 189219          1 817514
    ## 4 a       t       -+-        0.199   39112          1 870176
    ## 5 c       g       --+        0.730   80434          1 928622
    ## 6 t       c       --+        0.913  185512          1 996120

## Step 2: Map the SNP location to GRCh37

For MAGMA gene analysis, the required population-specific reference data
is based on the GRCh37 genome build, all SNP location should be
converted to GRCh37 (Hg19). The example data file shown above is a GWAS
summary statistics file based on the GRCh38 genome build. To map the SNP
location to GRCh37, we can use the [LiftOver
tool](https://genome.ucsc.edu/cgi-bin/hgLiftOver) provided by UCSC. Here
is an example of how to convert the SNP location with the chain object
“hg38ToHg19.over.chain.gz” downloaded from the [UCSC
website](https://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/),
together with the rtrackler package in R.

``` r
#load the chain object
chainObject <- import.chain("data/ref/hg38ToHg19.over.chain") #your downloaded chain object

#Granges object with the bed information. The seqnames are named with "chr" prefix.
grObject <- GRanges(seqnames=paste0("chr", example_gwas$Chromosome),
                   ranges=IRanges(start=as.numeric(example_gwas$BP), end=as.numeric(example_gwas$BP)+1))

#lift over
hg19_res <- liftOver(grObject, chainObject)

#extract the result data frame, select the necessary columns
hg19_res_df <- as.data.frame(hg19_res) %>% #get the result data frame
  mutate(seqnames = gsub(seqnames, pattern = "chr",replacement = "")) %>% #change pattern back
  select(group, start, seqnames) %>%
  rename("start" = "BP", "seqnames" = "Chromosome") #rename columns

#get the 1 to 1 mapping
hg19_res_df <- hg19_res_df %>% group_by(group) %>% 
  filter(n()==1) %>%
  ungroup() 

#merge with the original data frame (the group column indicates the original position in the table)
example_gwas <- example_gwas %>%
  mutate(group = 1:n()) %>%  #original order 
  select(-Chromosome, -BP) %>% #remove the original column
  left_join(hg19_res_df , by = "group") %>% #left join annotation
  filter(!is.na(BP)) #filter unmapped data

#write out for later
write.table(example_gwas %>% mutate(BP = as.character(BP)),  #avoid numerical output 
            file = "data/gwas/sample_gwas/sample_gwas_hg19.txt", 
            sep="\t", col.names = T, row.names = F, quote = F)

head(example_gwas)
```

    ## # A tibble: 6 × 7
    ##   Allele1 Allele2 Direction P.value group     BP Chromosome
    ##   <chr>   <chr>   <chr>       <dbl> <int>  <int> <chr>     
    ## 1 t       c       ---        0.0541     1 731718 1         
    ## 2 a       g       +++        0.0315     2 752721 1         
    ## 3 t       c       +++        0.0159     3 752894 1         
    ## 4 a       t       -+-        0.199      4 805556 1         
    ## 5 c       g       --+        0.730      5 864002 1         
    ## 6 t       c       --+        0.913      6 931500 1

## Step 3: Map the SNP location to SNP ID

If the SNP ID is not included in the GWAS summary statistics file, we
can use the SNP synonyms file to map the SNP location to SNP ID. To
simplify the procedure, we implement a bash script in
`tools/annotate_SNP.sh` for the step, where
[bedtools](https://bedtools.readthedocs.io/en/latest/) is required. This
script will require the dbSNP file and the GWAS summary statistics file
as input, and output a new file with SNP ID included in the last column.
To view all specific parameters of the script, you can run the following
command in the terminal:

``` bash
bash src/tools/annotate_SNP.sh -h
```

    ## 
    ##     This script help you locate SNP and transform the locus (by chromosome id and location) to the corresponding variant rsid.
    ##     Only intersected SNPs will be printed out.
    ##     Note!!!! 1. All file path should be the same as the original files. 2. bedtools is required, so make sure it's in your environment or path
    ##     Here are the required parameters:
    ##     -g | GWAS summary statistics files. Or any other files that you would like to annotate SNPs to.
    ##     -o | Outpur file header. If the value is given as 'dir1/dir2/output', then the final output file will be 'dir1/dir2/output.annot.table.txt'. If it is not specified, the output file will be the same as the input file with the suffix ignored.
    ##     -c | Chromosome column, i.e. the column with with chromosome information be (in the table). The value should be specified or by default it will be 1.
    ##     -p | Position column, i.e. The column with the the position information (in the table), The value should be specified or by default it will be 2.
    ##     -b | Basis of chromosome numbering, to specify the if first position on a chromosome is labeled as '1'(1-base), or '0'(0-base). The default parameters are start with 1. #only these two are available!
    ##     -k | Keep all unmapped SNPs which do not exist in the dbSNP file or not. By default it is True. (You can choose between 'True' or 'False')
    ##     
    ##     optional parameters:
    ##     -v | VCF file for the dbsnp version you used. For default it will be located in 'data/ref/Hg19.dbsnp151.vcf'. This can be downloaded online.
    ##     -n | Number of tasks to split the job. The default value is 1. If the memory is not enough, you can split the job into several parts but it will take more time.
    ##     -t | Customized bedtools command. By default it is 'bedtools'. If you do not have it in your environment, you can specify the path to the bedtools command.
    ##     -w | The working directory. By default it is the current directory. If you want to specify the working directory, you can use this parameter.
    ##     -h | Help message.
    ## 

Here is an example of how to run the script with the example data file:

``` bash
bash src/tools/annotate_SNP.sh -g data/gwas/sample_gwas/sample_gwas_hg19.txt -o data/gwas/sample_gwas/sample_gwas_hg19 \
-c 6 -p 7 -v data/ref/Hg19.dbsnp151.vcf -t bin/bedtools
```

Note that the script may takes hours to run and consumes more than 100
GB of memory. To alleviate memory usage (but increase the running time),
the `-n` option may be used to sequentially split the task into several.
The annotated table will be printed with the path and header specified
by the `-o` option.

``` bash
head data/gwas/sample_gwas/sample_gwas_hg19.annot.table
```

    ## Allele1  Allele2 Direction   P.value group   Chromosome  BP  marker.name
    ## t    c   --- 0.0541  1   1   731718 rs58276399
    ## a    g   +++ 0.03149 2   1   752721 rs3131972
    ## t    c   +++ 0.01587 3   1   752894 rs3131971
    ## a    t   -+- 0.1987  4   1   805556 rs72631880
    ## c    g   --+ 0.7303  5   1   864002 rs1806501
    ## t    c   --+ 0.913   6   1   931500 rs74718486
    ## a    g   --+ 0.4902  7   1   938125 rs2710868
    ## a    g   --+ 0.785   8   1   941334 rs57683598
    ## t    c   +++ 0.1863  9   1   948921 rs15842

## Step 4: Generate the input file for MAGMA

The software takes two input files generated from the previously
processed (or unprocessed but with the required columns) GWAS summary
statistics file: a SNP location file and a P value file. You can simply
use the `awk` command to extract the information, for example:

``` bash
#print rsID, chromosome, BP but not colum names:
#starting from the second line (NR>1), field 8, 6, 7
awk '{if (NR>1) print $8, $6, $7}' data/gwas/sample_gwas/sample_gwas_hg19.annot.table > data/gwas/sample_gwas/sample_gwas_hg19.snp_loc.txt
#print rsID, Pvalue but not colum names:
#starting from the second line (NR>1), field 8, 4
awk '{if (NR>1) print $8, $4}' data/gwas/sample_gwas/sample_gwas_hg19.annot.table > data/gwas/sample_gwas/sample_gwas_hg19.p_val.txt
```

or in R:

``` r
anno_tbl <- read.table("data/gwas/sample_gwas/sample_gwas_hg19.annot.table", header=T)
anno_tbl %>% select(marker.name, Chromosome, BP) %>% write.table(file = "data/gwas/sample_gwas/sample_gwas_hg19.snp_loc.txt", col.names = F, row.names = F, quote = F, sep="\t")
anno_tbl %>% select(marker.name, P.value) %>% write.table(file = "data/gwas/sample_gwas/sample_gwas_hg19.p_val.txt", col.names = F, row.names = F, quote = F, sep="\t")
```

Data preview:

``` bash
head  data/gwas/sample_gwas/sample_gwas_hg19.snp_loc.txt
```

    ## rs58276399   1   731718
    ## rs3131972    1   752721
    ## rs3131971    1   752894
    ## rs72631880   1   805556
    ## rs1806501    1   864002
    ## rs74718486   1   931500
    ## rs2710868    1   938125
    ## rs57683598   1   941334
    ## rs15842  1   948921
    ## rs9442364    1   970215

## Step 5: Run MAGMA software

There are two steps: gene annotation and gene analysis. You can use our
implemented bash script `tools/magma_gene_zscore_analysis.sh` to run the
two steps together. Another paramter to specificy is the gene analysis
window size, which is set as 35 kilobase (kb) upstream to 10 kb
downstream by default, as recommend by [the previous
study](https://www.nature.com/articles/s41588-020-0610-9). Although we
have shown that output files from various window sizes make only a minor
difference in the results of *seismic* analysis, you can still change
the window size by specifying the parameter `--window-size` in the
script. To view all specific parameters of the script, you can run the
following command in the terminal:

``` bash
bash src/tools/magma_gene_zscore_analysis.sh -h
```

    ## 
    ## Usage: src/tools/magma_gene_zscore_analysis.sh -l SNP_LOC_FILE -p SNP_P_VALUE_FILE -o OUTPUT_DIR -s COHORT_SIZE -m MAGMA_PATH -g GENE_LOC_FILE -b B_FILE [-w WINDOW_SIZES] [-h]
    ## 
    ## This script performs MAGMA annotation and analysis on SNP data.
    ## 
    ## Required parameters:
    ##     -l | SNP location file
    ##     -p | SNP p-value file
    ##     -o | Output directory
    ##     -s | Cohort size (use 'cX' format to specify column X containing sample size)
    ##     -m | Path to MAGMA binary
    ##     -g | Path to gene location file
    ##     -b | Path to g1000_eur files (without extension)
    ## 
    ## Optional parameters:
    ##     -w | Window sizes (default: 35,10). Use colon-separated list for multiple sizes, e.g., '35,10:15,20:40,10'
    ##     -h | Display this help message

``` bash
bash src/tools/magma_gene_zscore_analysis.sh -l data/gwas/sample_gwas/sample_gwas_hg19.snp_loc.txt \
-p data/gwas/sample_gwas/sample_gwas_hg19.p_val.txt \
-o data/gwas/sample_gwas \
-s 100000 -m bin/magma/magma -g data/ref/magma/NCBI37.3/NCBI37.3.gene.loc -b data/ref/magma/g1000_eur/g1000_eur 
```

Besides specifying these reference file paths and the curated files from
the summary statistics, the total size of the cohort should also be
specified. Note that sometimes the cohort size for each SNP may be
different, for example, due to some quality control procedure. In that
case, a separate column in the p-value file should be added and the
input `-s` argument should tell the script the order of the column (for
example, the third column) describing the cohort size. By specifying
these arguments, the gene-level MAGMA z-score will be printed out in the
output directory with header similar to the snp location file. Here the
output file is `data/gwas/sample_gwas/sample_gwas_hg19.35.10.genes.out`.
In addition, please do remember to check the MAGMA running log.
Sometimes when the data are not in the correct format, there will be
very few genes and SNPs annotated (for example, less than 10% of the
SNPs).
