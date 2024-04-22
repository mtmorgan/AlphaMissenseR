#' Default ClinVar Data
#'
#' Derived from the supplemental table of the AlphaMissense
#'      [\[2023\]](https://www.science.org/doi/10.1126/science.adg7492) paper.
#'
#' @format dataframe with 82872 rows and 5 variables:
#' \describe{
#'      \item{cv_variant_id}{ClinVar variant identifer.}
#'      \item{uniprot_id}{UniProt accession identifier.}
#'      \item{transcript_id}{Ensembl transcript identifier.}
#'      \item{protein_variant}{Protein variant identifier.}
#'      \item{cv_class}{Binary ClinVar classification.}
#'      }
#'
#' @examples
#' data(clinvar_data)
"clinvar_data"
