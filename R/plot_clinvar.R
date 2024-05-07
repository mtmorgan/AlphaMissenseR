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
#' For `am_table`, columns must include:
#'
#' - `uniprot_id`: UniProt accession identifier.
#' - `protein_variant`: variant identifier string, with protein
#'  position in the middle, and the reference and mutant amino acid
#'  residues to the left and right of the position, respectively.
#' - `am_class`: AlphaMissense classification of either "benign",
#'  "ambiguous", or "pathogenic".
#' - `am_pathogenicity`: AlphaMissense predicted score.
#'
#' For `cv_table`, columns must include:
#'
#'- `uniprot_id`: UniProt accession identifier, matching `am_table`.
#' - `protein_variant`: variant identifier string, matching `am_table`
#' format.
#' - `cv_class`: binary ClinVar classification of 0 (benign) or 1
#' (pathogenic).
#'
#' @param uniprotId a string with a valid UniProt accession
#'     identifier.
#'
#' @param am_table a table derived from `am_data("aa_substitution")`.
#'     Alternatively, a user-defined tibble or data.frame.
#'
#' @param cv_table a table containing ClinVar information derived from
#'     the supplemental table of the AlphaMissense
#'     [\[2023\]](https://www.science.org/doi/10.1126/science.adg7492)
#'     paper.  Alternatively, a user-defined tibble or data.frame.
#'
#' @return
#'
#' `plot_clinvar()` returns a `ggplot` object which overlays ClinVar
#' classifications onto AlphaMissense predicted scores. Blue, gray, and red
#' colors represent pathogenicity classifications for "likely benign",
#' "ambiguous", or "likely pathogenic", respectively. Large, bolded points
#' are ClinVar variants colored according to their clinical classification,
#' while smaller points in the background are AlphaMissense predictions.
#'
#' @examples
#'
#' data(clinvar_data)
#'
#' am_table <- db_connect() |>
#'             tbl("aa_substitutions") |>
#'             filter(uniprot_id == "P37023") |>
#'             dplyr::collect()
#'
#' plot_clinvar(uniprotId = "P37023",
#'              am_table = am_table,
#'              cv_table = clinvar_data)
#'
#' @importFrom BiocBaseUtils isCharacter
#' @import ggplot2
#' @import dplyr
#' @importFrom utils data
#'
#' @export
plot_clinvar <-
    function(uniprotId, am_table, cv_table)
{
    stopifnot(isCharacter(uniprotId))

    # am_table
    # If am_table not provided, load default
    if (missing(am_table)) {
        message("'am_table' not provided, using default ",
                "'am_data(\"aa_substitution\")' table accessed through ",
                "the AlphaMissenseR package.")

        am_table <- db_connect() |>
            tbl("aa_substitutions") |>
            filter(.data$uniprot_id == uniprotId) |>
            dplyr::collect()

    } else {
        # Take user-defined am_table and filter for the uniprotId
        am_table <- am_table |>
            filter(.data$uniprot_id == uniprotId)
    }

    # cv_table
    # If cv_table not provided, load default
    if (missing(cv_table)) {
        message("'cv_table' not provided, using default ",
                "ClinVar dataset in AlphaMissenseR package")

        data_env <- new.env(parent = emptyenv())
        data("clinvar_data", envir = data_env, package = "AlphaMissenseR")
        clinvar_data <- data_env[["clinvar_data"]]

        cv_table <-  clinvar_data |>
            filter(.data$uniprot_id == uniprotId)

    } else {
    # Take user-defined cv_table and filter for the uniprotId
        cv_table <- cv_table |>
            filter(.data$uniprot_id == uniprotId)
    }

    # Check that protein exists in ClinVar and AlphaMissense data
    if (!nrow(cv_table)) {
        stop(
            "No ClinVar information found for the protein accession. ",
            "Check that the UniProt ID is correct."
        )
    }

    if (!nrow(am_table)) {
        stop(
            "No AlphaMissense information found for the protein accession.",
            " Check that the UniProt ID is correct."
        )
    }

    # Validity check for am_table and cv_table
    am_required_columns <- c("uniprot_id", "protein_variant",
                            "am_class", "am_pathogenicity")
    stopifnot(
        inherits(am_table, "tbl") || inherits(am_table, "data.frame"),
        all(am_required_columns %in% colnames(am_table))
    )

    cv_required_columns <- c("uniprot_id", "protein_variant", "cv_class")
    stopifnot(
        inherits(cv_table, "tbl") || inherits(cv_table, "data.frame"),
        all(cv_required_columns %in% colnames(cv_table))
    )

    # Join databases by protein_variant
    am_table <- am_table |>
        mutate(
            aa_pos = as.numeric(
                gsub(".*?([0-9]+).*", "\\1", .data$protein_variant)
            )
        )

    c_combo <- left_join(
        am_table,
        cv_table,
        by = c('uniprot_id', 'protein_variant')
    )

    ##### PLOTTING #####

    # Grab the thresholds for AM pathogenicity to plot
    am_cutoff <- c_combo |>
        filter(.data$am_class == "ambiguous") |>
        select(.data$am_pathogenicity) |>
        summarise(
            path_cutoff = max(.data$am_pathogenicity),
            benign_cutoff = min(.data$am_pathogenicity)
        )

    # Add color code matching AM and CV labels
    c_combo <- c_combo |>
        mutate(code_color = case_when(
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
                )
            ) |>
        arrange(.data$code_color)

    # Plot sequence window
    plot <- c_combo |>
        ggplot(aes(.data$aa_pos, .data$am_pathogenicity)) +
                geom_point(
                aes(
                    shape = .data$code_color,
                    color = .data$code_color,
                    size = .data$code_color,
                    fill = .data$code_color,
                    stroke = .data$code_color
                )
            ) +
            scale_shape_manual(
                values = c(19, 19, 19, 21, 21)
            ) +
            scale_discrete_manual(
                aesthetics = "stroke",
                values = c(0, 0, 0, 1.5, 1.5)
            ) +
            scale_size_manual(
                values = c(2, 2, 2, 4, 4)
            ) +
            scale_color_manual(
                values = c("gray", "#89d5f5", "#f56c6c", "black", "black")
            ) +
            scale_fill_manual(
                values = c("gray", "#89d5f5", "#f56c6c", "#007cb0", "#c70606")
            ) +
            geom_hline(
                yintercept = am_cutoff$path_cutoff, linetype = 2,
                color = "#c70606"
            ) +
            geom_hline(
                yintercept = am_cutoff$benign_cutoff,
                linetype = 2, color = "#007cb0"
            ) +
            labs(title = paste0("UniProt ID: ", uniprotId)) +
            xlab("amino acid position") +
            ylab("AlphaMissense score") +
            theme_classic()

    plot +
        theme(
            axis.text.x = element_text(size = 16),
            axis.text.y = element_text(size = 16),
            axis.title.y = element_text(size = 16),
            axis.title.x = element_text(size = 16),
            legend.title = element_blank(),
            legend.text = element_text(size = 11)
        )
}
