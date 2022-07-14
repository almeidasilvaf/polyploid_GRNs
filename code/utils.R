

#' Create a gene expression matrix from EBI's Expression Atlas
#'
#' @param id Character of experiment ID (e.g., "E-GEOD-61857").
#'
#' @return A gene expression with genes in rows and samples in columns.
#' @importFrom readr read_tsv
#' @importFrom purrr reduce
#' @importFrom stringr str_count
#' @importFrom tidyr separate
#' @examples
#' id <- "E-GEOD-61857"
#' exp <- ebi2exp(id)
ebi2exp <- function(id = NULL) {

    url <- paste0(
        "http://ftp.ebi.ac.uk/pub/databases/microarray/data/atlas/experiments/",
        id, "/", id, "-tpms.tsv"
    )
    exp <- read.csv(url, comment.char = "#", sep = "\t")

    # First column to rownames, then delete columns 1 and 2
    rownames(exp) <- exp[,1]
    exp[, c(1,2)] <- NULL

    sep_exp <- purrr::reduce(
        seq_along(exp), .init = exp,
        ~ .x |>
            tidyr::separate(
                names(exp)[.y],
                sep = ",",
                into = paste0(names(exp)[.y], "_col_",
                              seq(1 + max(stringr::str_count(exp[[.y]], ',')))),
                fill = "right"
            )
    )
    # To numeric matrix
    sep_exp[] <- sapply(sep_exp, as.numeric)
    sep_exp <- as.matrix(sep_exp)

    # Remove non-coding RNAs
    sep_exp <- sep_exp[!startsWith(rownames(sep_exp), "ENSRNA"), ]
    return(sep_exp)
}


#' Get protein sequences for the longest isoform of a gene from Ensembl
#'
#' @param proteome An AAStringSet object.
#'
#' @return An AAStringSet with the longest isoform only.
#' @importFrom Biostrings width
ensembl_longest_isoform <- function(proteome = NULL) {

    pnames <- gsub(".*gene:", "", names(proteome))
    pnames <- gsub(" .*", "", pnames)

    names(proteome) <- pnames
    proteome <- proteome[order(Biostrings::width(proteome), decreasing = TRUE),]
    proteome <- proteome[!duplicated(names(proteome)), ]
    return(proteome)
}


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
        names(seq) <- gsub("\\.[0-9]$", "", names(seq))

        return(seq)
    })
    return(proteomes)
}


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
        names(seq) <- gsub("\\.[0-9]$", "", names(seq))

        return(seq)
    })
    return(cds)
}

#' Update node IDs of an edge list given a data frame with new IDs.
#'
#' @param edges A data frame representing the edge list.
#' @param id_conversion A 2-column data frame with old IDs in column 1 and
#' new IDs in column 2.
#'
#' @return The same edge list given as input, but with IDs changed.
#' @noRd
update_edges <- function(edges, id_conversion) {

    names(edges)[1:2] <- c("n1", "n2")
    names(id_conversion)[1:2] <- c("old_id", "node1")
    id_conversion <- id_conversion[, 1:2]

    # Columns to select
    sel_cols <- c("node1", "node2")
    if(ncol(edges) == 3) {
        third_col <- names(edges)[3]
        sel_cols <- c(sel_cols, third_col)
    }

    # Replace old names with new names
    new_edges <- dplyr::inner_join(
        edges, id_conversion, by = c("n1" = "old_id")
    ) %>%
        dplyr::inner_join(
            ., id_conversion %>% dplyr::rename(node2 = node1),
            by = c("n2" = "old_id")
        ) %>%
        dplyr::select(dplyr::all_of(sel_cols)) %>%
        as.data.frame()

    return(new_edges)
}


#' Convert protein IDs from STRINGdb to Ensembl gene IDs
#'
#' @param edgelist_ppi A data frame with protein 1 and protein 2 in columns 1
#' and 2, respectively. Additional columns are ignored.
#' @param alias A 3-column data frame with STRINGdb protein aliases in other
#' sources, which can be obtained in the Download section of the STRING
#' web interface. Columns 1, 2, and 3 represent STRING protein IDs, alias,
#' and source, respectively.
#' @param species Character scalar of species abbreviation with
#' first letter of genus and full specific epithet (e.g., "Athaliana").
#'
#' @return The same data frame passed as input, but with updated IDs.
#' @noRd
stringdb2ensembl <- function(edgelist_ppi, alias, species) {


    field_alias <- "BLAST_UniProt_GN_ORFNames"
    if(species == "Osativa") {
        field_alias <- "Ensembl_UniProt_GN"
    } else if(species == "Athaliana" | species == "Slycopersicum") {
        field_alias <- "Ensembl_gene"
    }
    # Get data frame of conversion between STRING IDs and Ensembl IDs
    names(alias) <- c("STRING", "alias", "source")
    names(edgelist_ppi)[1:2] <- c("protein1", "protein2")
    alias_ens <- alias[alias$source == field_alias, ]
    alias_ens <- alias_ens[!startsWith(alias_ens$alias, "LOC"), 1:2]

    if(species == "Vvinifera") {
        alias_ens <- readr::read_tsv(
            here::here("products", "result_files",
                       "vvinifera_mapping_string_ensembl.tsv"),
            show_col_types = FALSE
        )[, c(3,1)] %>% drop_na()
    }

    # Convert IDs
    edges <- update_edges(edgelist_ppi, alias_ens)

    # Handle special cases
    if(species == "Ptrichocarpa") {
        edges$node1 <- paste0(edges$node1, "v3")
        edges$node2 <- paste0(edges$node2, "v3")
    } else if(species == "Zmays") {
        edges$node1 <- gsub("ZEAMMB73_", "", edges$node1)
        edges$node2 <- gsub("ZEAMMB73_", "", edges$node2)
        v4_ids <- readr::read_tsv(
            "https://ftp.psb.ugent.be/pub/plaza/plaza_public_dicots_05/IdConversion/id_conversion.zma.csv.gz",
            show_col_types = FALSE, skip = 8
        ) %>% filter(id_type == "V4_identifier") %>% dplyr::select(3,1)
        edges <- update_edges(edges, as.data.frame(v4_ids))
    } else if(species == "Vvinifera") {
        ids <- readr::read_tsv(
            "https://ftp.psb.ugent.be/pub/plaza/plaza_public_dicots_05/IdConversion/id_conversion.vvi.csv.gz",
            show_col_types = FALSE, skip = 8
        )
    } else if(species == "Slycopersicum") {
        edges$node1 <- gsub("\\.[0-9]$", "", edges$node1)
        edges$node2 <- gsub("\\.[0-9]$", "", edges$node2)
    }

    return(edges)
}


