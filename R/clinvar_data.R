#' Default ClinVar Data
#'
#' Derived from the supplemental table of the AlphaMissense
#'      [\[2023\]](https://www.science.org/doi/10.1126/science.adg7492) paper.
#'
#' @format dataframe with 82872 rows and 5 variables:
#' \describe{
#'     \item{cv_variant_id}{ClinVar variant identifer.}
#'     \item{uniprot_id}{UniProt accession identifier.}
#'     \item{transcript_id}{Ensembl transcript identifier.}
#'     \item{protein_variant}{Protein variant identifier.}
#'     \item{cv_class}{Binary ClinVar class. 0 for benign or 1 for pathogenic.}
#' }
#'
#' @export
clinvar_data <-
    function(record = ALPHAMISSENSE_RECORD, bfc = BiocFileCache())
{
    db_rname <- paste0("AlphaMissense_", record)
    db_tbl_name <- "clinvar"
    if (!NROW(bfcquery(bfc, db_rname))) {
        spdl::info("creating AlphaMissense database for record '{}'", record)
        bfcnew(bfc, db_rname)
    }
    fpath <- system.file(
        package = "AlphaMissenseR", "extdata", "science.adg7492_data_s5.csv.gz"
    )
    am_data_import_csv(record, bfc, db_tbl_name, fpath = fpath, delim = ",")
}
