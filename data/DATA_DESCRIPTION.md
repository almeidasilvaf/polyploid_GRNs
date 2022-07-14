Data description
================

# Overview

In this project, we used the following species and associated RNA-seq
projects from EBI’s Expression Atlas:

| Species                    | Genes | Version                  | Source             | EA_Project   | N   |
|:---------------------------|------:|:-------------------------|:-------------------|:-------------|:----|
| Oryza sativa spp. japonica | 35775 | IRGSP-1.0                | Ensembl Plants v53 | E-MTAB-2037  | 35  |
| Zea mays                   | 39756 | Zm-B73-REFERENCE-NAM-5.0 | Ensembl Plants v53 | E-GEOD-50191 | 115 |
| Vitis vinifera             | 35134 | PN40024.v4               | Ensembl Plants v53 | E-GEOD-62744 | 100 |
| Glycine max                | 55897 | Glycine_max_v2.1         | Ensembl Plants v53 | E-GEOD-61857 | 40  |
| Solanum lycopersicum       | 34429 | SL3.0                    | Ensembl Plants v53 | E-MTAB-4813  | 90  |
| Populus trichocarpa        | 41335 | Pop_tri_v3               | Ensembl Plants v53 | E-GEOD-81077 | 40  |
| Arabidopsis thaliana       | 27628 | TAIR10                   | Ensembl Plants v53 | E-MTAB-7978  | 270 |

# Data files

Below, you can see a list of files in this directory, their
descriptions, and (when applicable) code to obtain it.

``` r
library(here)
```

## expression_data.rda

We will store gene expression matrices in a list for each species.
Expression values are represented in TPM and quantile-normalized.

``` r
#----Source function to get expression matrix from EBI--------------------------
source(here("code", "utils.R"))

#----Get expression matrices----------------------------------------------------
projects <- c(
    "E-MTAB-2037", "E-GEOD-50191", "E-GEOD-62744", "E-GEOD-61857",
    "E-MTAB-4813", "E-GEOD-81077", "E-MTAB-7978"
)

expression_data <- lapply(projects, ebi2exp)
names(expression_data) <- c(
    "Osativa", "Zmays", "Vvinifera", "Gmax", "Slycopersicum", 
    "Ptrichocarpa", "Athaliana"
)

#----Remove genes with median TPM < 1-------------------------------------------
expression_data <- lapply(expression_data, function(x) {
    exp <- BioNERO::remove_nonexp(x, min_exp = 1)
    return(exp)
})

#----Clean S. lycopersicum names------------------------------------------------
rownames(expression_data$Slycopersicum) <- gsub(
    "\\.[0-9]$", "", rownames(expression_data$Slycopersicum)
)

# Save data
save(
    expression_data,
    file = here("data", "expression_data.rda"),
    compress = "xz"
)
```

## proteomes.rda

Proteomes for each species stored as a list of `AAStringSet` objects.
Here, we only used the translated sequences of primary transcripts. For
consistency with the expression data, data were downloaded from Ensembl
Plants.

``` r
#----Get urls to proteomes files------------------------------------------------
proteome_urls <- c(
    Osativa = "http://ftp.ebi.ac.uk/ensemblgenomes/pub/release-53/plants/fasta/oryza_sativa/pep/Oryza_sativa.IRGSP-1.0.pep.all.fa.gz",
    Zmays = "http://ftp.ebi.ac.uk/ensemblgenomes/pub/release-53/plants/fasta/zea_mays/pep/Zea_mays.Zm-B73-REFERENCE-NAM-5.0.pep.all.fa.gz",
    Vvinifera = "http://ftp.ebi.ac.uk/ensemblgenomes/pub/release-53/plants/fasta/vitis_vinifera/pep/Vitis_vinifera.PN40024.v4.pep.all.fa.gz",
    Gmax = "http://ftp.ebi.ac.uk/ensemblgenomes/pub/release-53/plants/fasta/glycine_max/pep/Glycine_max.Glycine_max_v2.1.pep.all.fa.gz",
    Slycopersicum = "http://ftp.ebi.ac.uk/ensemblgenomes/pub/release-53/plants/fasta/solanum_lycopersicum/pep/Solanum_lycopersicum.SL3.0.pep.all.fa.gz",
    Ptrichocarpa = "http://ftp.ebi.ac.uk/ensemblgenomes/pub/release-53/plants/fasta/populus_trichocarpa/pep/Populus_trichocarpa.Pop_tri_v3.pep.all.fa.gz",
    Athaliana = "http://ftp.ebi.ac.uk/ensemblgenomes/pub/release-53/plants/fasta/arabidopsis_thaliana/pep/Arabidopsis_thaliana.TAIR10.pep.all.fa.gz"
)

#----Read proteomes as AAStringSet objects--------------------------------------
proteomes <- lapply(proteome_urls, function(x) {
    seq <- Biostrings::readAAStringSet(x)
    
    # Keep only longest isoform for each gene
    seq <- ensembl_longest_isoform(seq)
    return(seq)
})

#----Save object----------------------------------------------------------------
save(
    proteomes,
    file = here("data", "proteomes.rda"),
    compress = "xz"
)
```

Here, to save storage, we will define a function to obtain the object
`proteomes` on the fly whenever we need to use it.

``` r
get_proteomes <- function() {
    proteome_urls <- c(
        Osativa = "http://ftp.ebi.ac.uk/ensemblgenomes/pub/release-53/plants/fasta/oryza_sativa/pep/Oryza_sativa.IRGSP-1.0.pep.all.fa.gz",
        Zmays = "http://ftp.ebi.ac.uk/ensemblgenomes/pub/release-53/plants/fasta/zea_mays/pep/Zea_mays.Zm-B73-REFERENCE-NAM-5.0.pep.all.fa.gz",
        Vvinifera = "http://ftp.ebi.ac.uk/ensemblgenomes/pub/release-53/plants/fasta/vitis_vinifera/pep/Vitis_vinifera.PN40024.v4.pep.all.fa.gz",
        Gmax = "http://ftp.ebi.ac.uk/ensemblgenomes/pub/release-53/plants/fasta/glycine_max/pep/Glycine_max.Glycine_max_v2.1.pep.all.fa.gz",
        Slycopersicum = "http://ftp.ebi.ac.uk/ensemblgenomes/pub/release-53/plants/fasta/solanum_lycopersicum/pep/Solanum_lycopersicum.SL3.0.pep.all.fa.gz",
        Ptrichocarpa = "http://ftp.ebi.ac.uk/ensemblgenomes/pub/release-53/plants/fasta/populus_trichocarpa/pep/Populus_trichocarpa.Pop_tri_v3.pep.all.fa.gz",
        Athaliana = "http://ftp.ebi.ac.uk/ensemblgenomes/pub/release-53/plants/fasta/arabidopsis_thaliana/pep/Arabidopsis_thaliana.TAIR10.pep.all.fa.gz"
    )
    
    proteomes <- lapply(proteome_urls, function(x) {
        seq <- Biostrings::readAAStringSet(x)
        
        # Keep only longest isoform for each gene
        seq <- ensembl_longest_isoform(seq)
        return(seq)
    })
    return(proteomes)
}
```

This function was inserted in *utils.R* for reuse.

## tfs.rda

This is a list of data frames with IDs of transcript factor
(TF)-encoding genes in the first column and TF family in the second
column.

``` r
#----Load proteomes-------------------------------------------------------------
load(here("data", "proteomes.rda"))

#----Identify and classify TFs--------------------------------------------------
library(tfhunter)
tfs <- lapply(proteomes, function(x) {
    tf_annot <- annotate_pfam(x, mode = "local")
    tf_class <- classify_tfs(tf_annot)
    return(tf_class)
})

#----Clean S. lycopersicum IDs--------------------------------------------------
tfs$Slycopersicum$Gene <- gsub(
    "\\.[0-9]$", "", tfs$Slycopersicum$Gene
)

#----Save object----------------------------------------------------------------
save(
    tfs,
    file = here("data", "tfs.rda"),
    compress = "xz"
)
```

## cds.rda

CDS for each species stored as a list of `DNAStringSet` objects. Here,
we only used the CDS of primary transcripts. For consistency with the
expression data, data were downloaded from Ensembl Plants.

``` r
#----Get urls to CDS files------------------------------------------------
cds_urls <- c(
    Osativa = "http://ftp.ebi.ac.uk/ensemblgenomes/pub/release-53/plants/fasta/oryza_sativa/cds/Oryza_sativa.IRGSP-1.0.cds.all.fa.gz",
    Zmays = "http://ftp.ebi.ac.uk/ensemblgenomes/pub/release-53/plants/fasta/zea_mays/cds/Zea_mays.Zm-B73-REFERENCE-NAM-5.0.cds.all.fa.gz",
    Vvinifera = "http://ftp.ebi.ac.uk/ensemblgenomes/pub/release-53/plants/fasta/vitis_vinifera/cds/Vitis_vinifera.PN40024.v4.cds.all.fa.gz",
    Gmax = "http://ftp.ebi.ac.uk/ensemblgenomes/pub/release-53/plants/fasta/glycine_max/cds/Glycine_max.Glycine_max_v2.1.cds.all.fa.gz",
    Slycopersicum = "http://ftp.ebi.ac.uk/ensemblgenomes/pub/release-53/plants/fasta/solanum_lycopersicum/cds/Solanum_lycopersicum.SL3.0.cds.all.fa.gz",
    Ptrichocarpa = "http://ftp.ebi.ac.uk/ensemblgenomes/pub/release-53/plants/fasta/populus_trichocarpa/cds/Populus_trichocarpa.Pop_tri_v3.cds.all.fa.gz",
    Athaliana = "http://ftp.ebi.ac.uk/ensemblgenomes/pub/release-53/plants/fasta/arabidopsis_thaliana/cds/Arabidopsis_thaliana.TAIR10.cds.all.fa.gz"
)

#----Read proteomes as AAStringSet objects--------------------------------------
cds <- lapply(cds_urls, function(x) {
    seq <- Biostrings::readDNAStringSet(x)
    
    # Keep only longest isoform for each gene
    seq <- ensembl_longest_isoform(seq)
    return(seq)
})

#----Save object----------------------------------------------------------------
save(
    cds,
    file = here("data", "cds.rda"),
    compress = "xz"
)
```

Here, to save storage, we will define a function to create the object
`cds` on the fly whenever we need to use it.

``` r
get_cds <- function() {
    cds_urls <- c(
        Osativa = "http://ftp.ebi.ac.uk/ensemblgenomes/pub/release-53/plants/fasta/oryza_sativa/cds/Oryza_sativa.IRGSP-1.0.cds.all.fa.gz",
        Zmays = "http://ftp.ebi.ac.uk/ensemblgenomes/pub/release-53/plants/fasta/zea_mays/cds/Zea_mays.Zm-B73-REFERENCE-NAM-5.0.cds.all.fa.gz",
        Vvinifera = "http://ftp.ebi.ac.uk/ensemblgenomes/pub/release-53/plants/fasta/vitis_vinifera/cds/Vitis_vinifera.PN40024.v4.cds.all.fa.gz",
        Gmax = "http://ftp.ebi.ac.uk/ensemblgenomes/pub/release-53/plants/fasta/glycine_max/cds/Glycine_max.Glycine_max_v2.1.cds.all.fa.gz",
        Slycopersicum = "http://ftp.ebi.ac.uk/ensemblgenomes/pub/release-53/plants/fasta/solanum_lycopersicum/cds/Solanum_lycopersicum.SL3.0.cds.all.fa.gz",
        Ptrichocarpa = "http://ftp.ebi.ac.uk/ensemblgenomes/pub/release-53/plants/fasta/populus_trichocarpa/cds/Populus_trichocarpa.Pop_tri_v3.cds.all.fa.gz",
        Athaliana = "http://ftp.ebi.ac.uk/ensemblgenomes/pub/release-53/plants/fasta/arabidopsis_thaliana/cds/Arabidopsis_thaliana.TAIR10.cds.all.fa.gz"
    )
    
    cds <- lapply(cds_urls, function(x) {
        seq <- Biostrings::readDNAStringSet(x)
        
        # Keep only longest isoform for each gene
        seq <- ensembl_longest_isoform(seq)
        return(seq)
    })
    return(cds)
}
```

This function was inserted in *utils.R* for reuse.

## annotation.rda

The gene annotation will be stored in a list of GRanges objects. Again,
for consistency, data will be downloaded from Ensembl Plants.

``` r
#----Get urls to CDS files------------------------------------------------
annot_urls <- c(
    Osativa = "http://ftp.ebi.ac.uk/ensemblgenomes/pub/release-53/plants/gtf/oryza_sativa/Oryza_sativa.IRGSP-1.0.53.gtf.gz",
    Zmays = "http://ftp.ebi.ac.uk/ensemblgenomes/pub/release-53/plants/gtf/zea_mays/Zea_mays.Zm-B73-REFERENCE-NAM-5.0.53.gtf.gz",
    Vvinifera = "http://ftp.ebi.ac.uk/ensemblgenomes/pub/release-53/plants/gtf/vitis_vinifera/Vitis_vinifera.PN40024.v4.53.gtf.gz",
    Gmax = "http://ftp.ebi.ac.uk/ensemblgenomes/pub/release-53/plants/gtf/glycine_max/Glycine_max.Glycine_max_v2.1.53.gtf.gz",
    Slycopersicum = "http://ftp.ebi.ac.uk/ensemblgenomes/pub/release-53/plants/gtf/solanum_lycopersicum/Solanum_lycopersicum.SL3.0.53.gtf.gz",
    Ptrichocarpa = "http://ftp.ebi.ac.uk/ensemblgenomes/pub/release-53/plants/gtf/populus_trichocarpa/Populus_trichocarpa.Pop_tri_v3.53.gtf.gz",
    Athaliana = "http://ftp.ebi.ac.uk/ensemblgenomes/pub/release-53/plants/gtf/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.53.gtf.gz"
)

#----Read proteomes as AAStringSet objects--------------------------------------
annotation <- lapply(annot_urls, function(x) {
    
    ranges <- rtracklayer::import(x)
    
    ranges <- ranges[ranges$gene_biotype == "protein_coding"]
    ranges <- ranges[ranges$type == "gene"]
    ranges$score <- NULL
    ranges$phase <- NULL
    ranges$gene_biotype <- NULL
    ranges$transcript_id <- NULL
    ranges$transcript_source <- NULL
    ranges$transcript_biotype <- NULL
    ranges$exon_id <- NULL
    ranges$exon_number <- NULL
    ranges$protein_id <- NULL
    ranges$transcript_name <- NULL
    ranges$gene_source <- NULL
    ranges$gene_id <- gsub("\\.[0-9]$", "", ranges$gene_id) 
    
    return(ranges)
})

#----Save object----------------------------------------------------------------
save(
    annotation,
    file = here("data", "annotation.rda"),
    compress = "xz"
)
```

## functional_annotation.rda

Functional annotation for whole gene set of each species will be
obtained from Ensembl with the Bioconductor package
*[biomaRt](https://bioconductor.org/packages/3.15/biomaRt)*.

``` r
load(here("data", "annotation.rda"))
library(biomaRt)

#----Choose Ensembl database----------------------------------------------------
ensembl <- useEnsemblGenomes(biomart = "plants_mart")

#----Choose a dataset-----------------------------------------------------------
datasets <- listDatasets(ensembl)

dataset_names <- list(
    Athaliana = "athaliana_eg_gene",
    Gmax = "gmax_eg_gene",
    Ptrichocarpa = "ptrichocarpa_eg_gene",
    Slycopersicum = "slycopersicum_eg_gene",
    Vvinifera = "vvinifera_eg_gene",
    Osativa = "osativa_eg_gene",
    Zmays = "zmays_eg_gene"
)

#----Get functional annotation--------------------------------------------------
functional_annotation <- lapply(dataset_names, function(x) {
    
    message("Working on dataset ", x)
    # Get Mart dataset
    ensembl_dataset <- useEnsemblGenomes(
        biomart = "plants_mart", 
        dataset = x
    )
    
    # GO
    gene_go <- getBM(
        attributes = c("ensembl_gene_id", "name_1006", "namespace_1003"),
        mart = ensembl_dataset
    )
    names(gene_go) <- c("gene", "GO", "type")
    gene_go <- gene_go[gene_go$GO != "", ]

    go_bp <- gene_go[gene_go$type == "biological_process", 1:2]
    go_mf <- gene_go[gene_go$type == "molecular_function", 1:2]
    go_cc <- gene_go[gene_go$type == "cellular_component", 1:2]

    # Interpro
    gene_interpro <- getBM(
        attributes = c("ensembl_gene_id", "interpro_description"),
        mart = ensembl_dataset
    )
    names(gene_interpro) <- c("gene", "interpro")
    gene_interpro <- gene_interpro[gene_interpro$interpro != "", ]
    

    # Combine results
    results <- list(
        GOBP = go_bp,
        GOMF = go_mf,
        GOCC = go_cc,
        InterPro = gene_interpro
    )
    if(x == "slycopersicum_eg_gene") {
        results <- lapply(results, function(x) {
            y <- x
            y$gene <- gsub("\\.[0-9]$", "", y$gene)
            return(y)
        })
    }
    return(results)
})

save(
    functional_annotation,
    file = here("data", "functional_annotation.rda"),
    compress = "xz"
)
```
