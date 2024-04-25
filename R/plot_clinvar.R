#' @rdname plot_clinvar
#'
#' @title Integrate ClinVar Labels with AlphaMissense Pathogenicity Scores
#'
#' @description `plot_clinvar()` integrates ClinVar classifications
#'      with AlphaMissense predicted scores derived from
#'      `am_data("aa_substitutions")` and returns a ggplot object for
#'      visualization.
#'
#' @param uniprotId a valid UniProt accession identifier
#'
#' @param am_table a table derived from `am_data("aa_substitution")`.
#'      Alternatively, a user-defined tibble or data.frame.
#'      Columns must include:
#'
#'- `uniprot_id`: UniProt accession identifier(s).
#'- `protein_variant`: variant identifier string, with the protein.
#'  position in the middle and the reference and mutant amino acid residues
#'  to the left and right of the position, respectively.
#' - `am_class`: AlphaMissense classification of either "benign",
#'  "ambiguous", or "pathogenic".
#' - `am_pathogenicity`: AlphaMissense predicted score.
#'
#' @param cv_table a table containing ClinVar information derived from the
#'      supplemental table of the AlphaMissense
#'      [\[2023\]](https://www.science.org/doi/10.1126/science.adg7492) paper.
#'      Alternatively, a user-defined tibble or data.frame.
#'      Columns must include:
#'
#'- `uniprot_id`: UniProt accession identifier(s), matching `am_table`.
#'- `protein_variant`: variant identifier string, matching `am_table` format.
#'- `cv_class`: binary ClinVar classification of 0 (benign) or 1 (pathogenic).
#'
#' @return `plot_clinvar()` returns a `ggplot` object which overlays
#'      ClinVar classifications onto AlphaMissense predicted scores for
#'      comparison and visualization.
#'
#' @examples
#' data(clinvar_data)
#' am_table <- db_connect() |>
#'             tbl("aa_substitutions") |>
#'             filter(uniprot_id == "P37023") |>
#'             dplyr::as_tibble()
#' plot_clinvar(uniprot_id = "P37023",
#'                 am_table = am_table,
#'                 cv_table = clinvar_data)
#'
#' @importFrom BiocBaseUtils isCharacter
#' @import ggplot2
#' @import dplyr
#'
#' @export
plot_clinvar <-
    function(uniprotId, am_table, cv_table)
    {
        # Validity checks
        if (isCharacter(uniprotId) == FALSE){
            stop(
                "uniprotId must be a string"
                )
        }


        # If am_table not provided, load default
        if (missing(am_table)){
            message("parameter 'am_table' not provided, using default ",
                    "'am_data(\"aa_substitution\")' table accessed through ",
                    "the AlphaMissenseR package.")

            am_table <- db_connect() |>
                tbl("aa_substitutions") |>
                filter(uniprot_id == uniprotId) |>
                dplyr::as_tibble()

        } else {
            # Take user-defined am_table and filter for the uniprotId
            am_table <- am_table |>
                filter(uniprot_id == uniprotId)
        }

        # If cv_table not provided, load default
        if (missing(cv_table)){
            message("parameter 'cv_table' not provided, using default ",
                    "ClinVar dataset in AlphaMissenseR package")

        data_env <- new.env(parent = emptyenv())
        data("clinvar_data", envir = data_env, package = "AlphaMissenseR")
        clinvar_data <- data_env[["clinvar_data"]]

        c_cv <-  clinvar_data |>
            filter(uniprot_id == uniprotId)

        } else {
        # Take user-defined cv_table and filter for the uniprotId
        c_cv <- cv_table |>
            filter(uniprot_id == uniprotId)
        }

        # Check that protein exists in ClinVar and AlphaMissense data
        if (nrow(c_cv) < 1){
            stop(
                "No ClinVar information found for the protein accession. ",
                "Check that the UniProt ID is correct."
            )
        }

        if (nrow(am_table) < 1) {
            stop(
                "No AlphaMissense information found for the protein accession.",
                " Check that the UniProt ID is correct."
            )
        }

        # Select relevant columns and join databases by protein_variant
        am_table <- am_table |>
            select(uniprot_id, protein_variant, am_pathogenicity, am_class)

        c_cv <- c_cv |>
            select(uniprot_id, protein_variant, cv_class)

        am_table <- am_table %>%
            mutate(aa_pos = as.numeric(gsub(".*?([0-9]+).*", "\\1",
                                                     protein_variant)))

        c_combo <- left_join(am_table, c_cv,
                             by = c('uniprot_id', 'protein_variant'))

        # Check if protein has multiple transcripts
        # unique_transcript_id <-
        #     c_combo |>
        #     filter(!is.na(transcript_id)) |>
        #     select(transcript_id) |>
        #     distinct() |>
        #     pull()
        #
        #  if (length(unique_transcript_id) > 1L){
        #      stop("Multiple transcripts detected")
        #  }
        #
        # # Mutate the new column with protein transcript
        # c_combo <- c_combo |>
        #     mutate(transcript_id = rep(unique_transcript_id))


        ##### PLOTTING #####

        # Grab the thresholds for AM pathogenicity to plot
        am_cutoff <- c_combo |>
            filter(am_class == "ambiguous") |>
            select(am_pathogenicity) |>
            summarise(
                path_cutoff = min(am_pathogenicity),
                benign_cutoff = max(am_pathogenicity)
            )

        # Add color code matching AM and CV labels
        c_combo <- c_combo |>
            mutate(code_color = case_when(
                !is.na(cv_class) & cv_class == "0" ~ "CV_benign",
                !is.na(cv_class) & cv_class == "1" ~ "CV_pathogenic",
                is.na(cv_class) & am_class == "pathogenic" ~ "AM_pathogenic",
                is.na(cv_class) & am_class == "benign" ~ "AM_benign",
                is.na(cv_class) & am_class == "ambiguous" ~ "AM_ambiguous"
                )
            )


        # Plot sequence window
        plot <- ggplot(c_combo |> arrange(code_color),
                       aes(aa_pos, am_pathogenicity)) +
                geom_point(
                aes(
                    shape = code_color,
                    color = code_color,
                    size = code_color,
                    fill = code_color,
                    stroke = code_color
                )
            ) +
            scale_shape_manual(values = c(19, 19, 19, 21, 21)) +
            scale_discrete_manual(aesthetics = "stroke",
                                  values = c(0, 0, 0, 1.5, 1.5)) +
            scale_size_manual(values = c(2, 2, 2, 4, 4)) +
            scale_color_manual(
                values = c("gray", "#89d5f5", "#f56c6c", "black", "black")
            ) +
            scale_fill_manual(
                values = c("gray", "#89d5f5", "#f56c6c", "#007cb0", "#c70606")
            ) +
            geom_hline(yintercept = am_cutoff$path_cutoff, linetype = 2,) +
            geom_hline(yintercept = am_cutoff$benign_cutoff,
                       linetype = 2) +
            labs(title = paste0("UniProt ID: ", uniprotId)) +
            xlab("amino acid position") +
            ylab("AlphaMissense score") +
            theme_classic()


        plot <- plot + theme(
            axis.text.x = element_text(size = 16),
            axis.text.y = element_text(size = 16),
            axis.title.y = element_text(size = 16),
            axis.title.x = element_text(size = 16),
            legend.position = "none"
        )

    plot

}
