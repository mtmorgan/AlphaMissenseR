#' @rdname plot_clinvar
#'
#' @title Integrate ClinVar Labels with AlphaMissense Pathogenicity Scores
#'
#' @description `plot_clinvar()` integrates ClinVar classifications with
#'    AlphaMissense predicted scores derived from `am_data("aa_substitutions")`
#'    and returns a ggplot object for visualization.
#'
#' @details
#'
#' For `alphamissense_table`, columns must include:
#'
#' - `uniprot_id`: UniProt accession identifier.
#' - `protein_variant`: variant identifier string, with protein
#'    position in the middle, and the reference and mutant amino acid
#'    residues to the left and right of the position, respectively.
#' - `am_class`: AlphaMissense classification of either "benign",
#'    "ambiguous", or "pathogenic".
#' - `am_pathogenicity`: AlphaMissense predicted score.
#'
#' For `clinvar_table`, columns must include:
#'
#' - `uniprot_id`: UniProt accession identifier, matching `alphamissense_table`.
#' - `protein_variant`: variant identifier string, matching
#'    `alphamissense_table` format.
#' - `cv_class`: binary ClinVar classification of 0 (benign) or 1 (pathogenic).
#'
#' @param uniprotId a string with a valid UniProt accession identifier.
#'
#' @param alphamissense_table a table containing AlphaMissense predictions for
#'    protein variants. By default, the table is derived from
#'    `am_data("aa_substitution")`. Alternatively, a user-defined
#'    [tibble::tibble()] or [base::data.frame()] can be supplied.
#'
#' @param clinvar_table a table containing ClinVar information. By default, the
#'    table is derived from the supplemental data of the AlphaMissense paper.
#'    Alternatively, a user-defined [tibble::tibble()] or [base::data.frame()]
#'    can be supplied.
#'
#' @return `plot_clinvar()` returns a `ggplot` object which overlays ClinVar
#'    classifications onto AlphaMissense predicted scores. Blue, gray, and red
#'    colors represent pathogenicity classifications for "likely benign",
#'    "ambiguous", or "likely pathogenic", respectively. Large, bolded points
#'    are ClinVar variants colored according to their clinical classification,
#'    while smaller points in the background are AlphaMissense predictions.
#'
#' @examples
#' data(clinvar_data)
#'
#' alphamissense_table <- am_data("aa_substitutions") |>
#'                            filter(uniprot_id == "P37023") |>
#'                            dplyr::collect()
#'
#' plot_clinvar(uniprotId = "P37023",
#'    alphamissense_table = alphamissense_table,
#'    clinvar_table = clinvar_data)
#'
#' @references Cheng et al.,
#' Accurate proteome-wide missense variant effect prediction with AlphaMissense.
#' \emph{Science} 381, eadg7492. DOI:10.1126/science.adg7492.
#'
#' @importFrom BiocBaseUtils isCharacter
#' @import ggplot2
#' @import dplyr
#' @importFrom utils data
#'
#' @export
plot_clinvar <-
    function(uniprotId, alphamissense_table, clinvar_table)
{
    stopifnot(isCharacter(uniprotId))

    # Load and check AM and CV tables
    alphamissense_table <- .check_am_table(am_table = alphamissense_table,
                                            uID = uniprotId)
    clinvar_table <- .check_cv_table(cv_table = clinvar_table,
                                    uID = uniprotId)

    # Validity check for alphamissense_table and clinvar_table
    if (!nrow(clinvar_table)) {
        stop(
            "No ClinVar information found for the protein accession. ",
            "Check that the UniProt ID is correct."
        )
    }

    if (!nrow(alphamissense_table)) {
        stop(
            "No AlphaMissense information found for the protein accession.",
            " Check that the UniProt ID is correct."
        )
    }

    am_required_columns <- c("uniprot_id", "protein_variant",
                            "am_class", "am_pathogenicity")
    stopifnot(
        inherits(alphamissense_table, "tbl") ||
            inherits(alphamissense_table, "data.frame"),
        all(am_required_columns %in% colnames(alphamissense_table))
    )

    cv_required_columns <- c("uniprot_id", "protein_variant", "cv_class")
    stopifnot(
        inherits(clinvar_table, "tbl") || inherits(clinvar_table, "data.frame"),
        all(cv_required_columns %in% colnames(clinvar_table))
    )

    # Join and process tables for plotting
    c_combo <- .process_data_for_plot_clinvar(am_table = alphamissense_table,
                                                cv_table = clinvar_table)

    # Grab the thresholds for AM pathogenicity to plot
    am_cutoff <- c_combo |>
        filter(.data$am_class == "ambiguous") |>
        select(.data$am_pathogenicity) |>
        summarise(
            path_cutoff = max(.data$am_pathogenicity),
            benign_cutoff = min(.data$am_pathogenicity)
        )

    # Plot sequence window
    .make_cv_plot(combo_table = c_combo, cutoff = am_cutoff, uId = uniprotId)
}
