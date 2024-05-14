#' @noRd
#'
#' Filter the AlphaMissense table with uniprotID
filter_am_table <-
    function(am_table, uID)
{
    ## Check if am_table is missing
    if (missing(am_table)) {
        message(
            "'alphamissense_table' not provided, using default ",
            "'am_data(\"aa_substitution\")' table accessed through ",
            "the AlphaMissenseR package."
        )
        am_table <- am_data("aa_substitutions")
    }

    ## Take alphamissense_table and filter for the uniprotId
    alphamissense_table <- am_table |>
            filter(.data$uniprot_id == uID) |>
            dplyr::collect()

    ## Check if table is empty after filtering
    ## This will work for a tibble or a data.frame
    if (!nrow(alphamissense_table)) {
        stop("No AlphaMissense information found for the protein accession.",
             " Check that the UniProt ID is correct.")
    }

    alphamissense_table
}


#' @noRd
#'
#' Filter the clinvar table with uniprot ID
filter_cv_table <-
    function(cv_table, uID)
{
    if (missing(cv_table)) {
        message("'clinvar_table' not provided, using default ",
                "ClinVar dataset in AlphaMissenseR package")

        data_env <- new.env(parent = emptyenv())
        data("clinvar_data", envir = data_env, package = "AlphaMissenseR")
        clinvar_data <- data_env[["clinvar_data"]]
    }

    ## Take clinvar_table and filter for the uniprotId
    clinvar_table <- cv_table |>
        filter(.data$uniprot_id == uID)

    ## Check if the table is empty after filtering
    if (!nrow(clinvar_table)) {
        stop("No ClinVar information found for the protein accession. ",
             "Check that the UniProt ID is correct.")
    }

    clinvar_table
}

#' @noRd
#'
#' Prepare data for the function plot_clinvar
prepare_data_for_plot_clinvar <-
    function(am_table, cv_table)
{
    ## grab amino acid positions
    am_table <- am_table |>
        mutate(
            aa_pos = as.numeric(
                gsub(".*?([0-9]+).*", "\\1", .data$protein_variant)
            )
        )

    ## join datasets
    combined_data <- left_join(
        am_table,
        cv_table,
        by = c('uniprot_id', 'protein_variant')
    )

    ## add color code matching AM and CV labels
    combined_data <-
        combined_data |>
        mutate(
            code_color = case_when(
                !is.na(.data$cv_class) & .data$cv_class == "0" ~
                    "CV benign",
                !is.na(.data$cv_class) & .data$cv_class == "1" ~
                    "CV pathogenic",
                is.na(.data$cv_class) & .data$am_class == "pathogenic" ~
                    "AM pathogenic",
                is.na(.data$cv_class) & .data$am_class == "benign" ~
                    "AM benign",
                is.na(.data$cv_class) & .data$am_class == "ambiguous" ~
                    "AM ambiguous"
            )) |>
        mutate_at(vars(code_color), factor) |>
        arrange(.data$code_color)

    ## Grab the thresholds for AM pathogenicity to plot
    combined_data |>
        group_by(.data$am_class) |>
        mutate(max = max(.data$am_pathogenicity, na.rm=TRUE),
               min = min(.data$am_pathogenicity, na.rm=TRUE))
}

#' @noRd
#'
#' Create a ClinVar plotting function using ggplot
create_clinvar_plot <-
    function(combined_table, uId)
{
    ## Create named vectors for all scale layers
    colScale <- scale_colour_manual(
        name = "code_color",
        values = c("AM ambiguous" = "gray",
                   "AM benign" = "#89d5f5",
                   "AM pathogenic" = "#f56c6c",
                   "CV benign" = "black",
                   "CV pathogenic" = "black"))

    fillScale <- scale_fill_manual(
        name = "code_color",
        values = c("AM ambiguous" = "gray",
                   "AM benign" = "#89d5f5",
                   "AM pathogenic" = "#f56c6c",
                   "CV benign" = "#007cb0",
                   "CV pathogenic" = "#c70606"))

    shapeScale <- scale_shape_manual(
        name = "code_color",
        values = c("AM ambiguous" = 19,
                   "AM benign" = 19,
                   "AM pathogenic" = 19,
                   "CV benign" = 21,
                   "CV pathogenic" = 21))

    sizeScale <- scale_size_manual(
        name = "code_color",
        values = c("AM ambiguous" = 2,
                   "AM benign" = 2,
                   "AM pathogenic" = 2,
                   "CV benign" = 4,
                   "CV pathogenic" = 4))

    strokeScale <- scale_discrete_manual(
        name = "code_color",
        aesthetics = "stroke",
        values = c("AM ambiguous" = 0,
                   "AM benign" = 0,
                   "AM pathogenic" = 0,
                   "CV benign" = 1.5,
                   "CV pathogenic" = 1.5))

    cv_plot <- combined_table |>
        ggplot(aes(.data$aa_pos, .data$am_pathogenicity)) +
        geom_point(
            aes(shape = .data$code_color,
                color = .data$code_color,
                size = .data$code_color,
                fill = .data$code_color,
                stroke = .data$code_color)
        ) +
        strokeScale +
        shapeScale +
        sizeScale +
        colScale +
        fillScale +
        geom_hline(
            yintercept = combined_table |>
                filter(am_class == "pathogenic") |>
                pull(min), linetype = 2,
            color = "#c70606"
        ) +
        geom_hline(
            yintercept = combined_table |>
                filter(am_class == "benign") |>
                pull(max),
            linetype = 2, color = "#007cb0"
        ) +
        labs(title = paste0("UniProt ID: ", uId)) +
        xlab("amino acid position") +
        ylab("AlphaMissense score") +
        theme_classic() +
        theme(
            axis.text.x = element_text(size = 16),
            axis.text.y = element_text(size = 16),
            axis.title.y = element_text(size = 16),
            axis.title.x = element_text(size = 16),
            legend.title = element_blank(),
            legend.text = element_text(size = 11)
        )

    cv_plot
}


#' @rdname plot_clinvar
#'
#' @title Integrate ClinVar Labels with AlphaMissense Pathogenicity
#'     Scores
#'
#' @description `plot_clinvar()` integrates ClinVar classifications
#'     with AlphaMissense predicted scores derived from
#'     `am_data("aa_substitutions")` and returns a ggplot object for
#'     visualization.
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
#' @param uniprotId a string with a valid UniProt accession
#'     identifier.
#'
#' @param alphamissense_table a table containing AlphaMissense
#'     predictions for protein variants. By default, the table is
#'     derived from `am_data("aa_substitution")`. Alternatively, a
#'     user-defined \code{\link{tibble}} or \code{\link{data.frame}}
#'     can be supplied.
#'
#' @param clinvar_table a table containing ClinVar information. By
#'     default, the table is derived from the supplemental data of the
#'     AlphaMissense paper.  Alternatively, a user-defined
#'     \code{\link{tibble}} or \code{\link{data.frame}} can be
#'     supplied.
#'
#' @return `plot_clinvar()` returns a `ggplot` object which overlays
#'     ClinVar classifications onto AlphaMissense predicted
#'     scores. Blue, gray, and red colors represent pathogenicity
#'     classifications for "likely benign", "ambiguous", or
#'     "likely pathogenic", respectively. Large, bolded points are
#'     ClinVar variants colored according to their clinical
#'     classification, while smaller points in the background are
#'     AlphaMissense predictions.
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
    ## Validate arguments
    stopifnot(isCharacter(uniprotId))

    ## Filter AM and CV tables with uniProtID
    alphamissense_table <- filter_am_table(
        am_table = alphamissense_table,
        uID = uniprotId
    )

    clinvar_table <- filter_cv_table(
        cv_table = clinvar_table,
        uID = uniprotId
    )

    ## Validate AM tables
    am_required_columns <- c("uniprot_id", "protein_variant",
                             "am_class", "am_pathogenicity")
    stopifnot(
        inherits(alphamissense_table, "tbl") ||
            inherits(alphamissense_table, "data.frame"),
        all(am_required_columns %in% colnames(alphamissense_table))
    )

    ## Validate CV table
    cv_required_columns <- c("uniprot_id", "protein_variant", "cv_class")
    stopifnot(
        inherits(clinvar_table, "tbl") ||
            inherits(clinvar_table, "data.frame"),
        all(cv_required_columns %in% colnames(clinvar_table))
    )

    ## Process tables for plotting
    combined_table <- prepare_data_for_plot_clinvar(
        am_table = alphamissense_table,
        cv_table = clinvar_table
    )

    ## Plot sequence window
    create_clinvar_plot(combined_table = combined_table, uId = uniprotId)
}
