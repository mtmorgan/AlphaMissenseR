#' @rdname plot_clinvar
#'
#' @title Integrate ClinVar Classifications with AlphaMissense Pathogenicity Scores
#'
#' @description `plot_clinvar()` integrates ClinVar classifications
#'      with AlphaMissense predicted scores derived from
#'      `am_data("aa_substitutions")` and returns a ggplot object for
#'      visualization.
#'
#' @param pID a valid UniProt accession identifier. # check Martin's code
#'
#' @param am_table a data.frame #link to description - derived from `am_data("aa_substitution")`.
#'      Alternatively, a user-defined tibble or dataframe.
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
#' @param cv_table a tibble or dataframe containing ClinVar information. By
#'      default, ClinVar information is derived from the supplemental table of
#'      the AlphaMissense
#'      [\[2023\]](https://www.science.org/doi/10.1126/science.adg7492) paper.
#'      Alternatively, a user-defined tibble or dataframe.
#'      Columns must include:
#'
#'- `uniprot_id`: UniProt accession identifier(s), matching AlphaMissense table.
#'- `protein_variant`: protein variant identifier string, matching AlphaMissense
#'  format.
#' - `cv_class`: binary values 0 (benign) and 1 (pathogenic) in ClinVar.
#'
#' @return `integrate_clinvar()` returns a `ggplot` object which overlays
#'      ClinVar pathogenicity annotations onto AlphaMissense predicted scores
#'      for comparison and visualization.
#'
#' @examples
#' data(clinvar_data)
#'
#' # get AM table with that pID first
#'
#' plot_clinvar(protein = "P37023", cv_table = clinvar_data, am_table = X)
#'
#' @importFrom BiocBaseUtils isCharacter
#' @import ggplot2
#' @import dplyr
#'
#' @export
integrate_clinvar <-
    function(protein, am_table, cv_table = clinvar_data)
    {
        # Validity checks
        stopifnot(isCharacter(protein))

        ### Need to add check on cv_table and am_table for required columns in
        ### user-defined scenario
        ### Check if user loaded

        # Load in default, inhouse ClinVar data
        data(clinvar_data)

        # Default tables
        c_cv <- clinvar_data |>
            filter(uniprot_id == protein)

        # Check if protein found in ClinVar
        if (nrow(c_cv) < 1){
            stop(
                "No ClinVar information found for the protein accession. ",
                "Check that the UniProt ID is correct."
            )
        }

        # Grab AM proteins that match CV
        c_am <- db_connect() |>
            tbl("aa_substitutions") |>
            filter(uniprot_id == protein) |>
            dplyr::as_tibble()

        # Check if protein found in AM
        if (nrow(c_am) < 1) {
            stop(
                "No AlphaMissense information found for the protein accession.",
                "Check that the UniProt ID is correct."
            )
        }

        # Select relevant columns from both data
        c_am <- c_am |>
            select(uniprot_id, protein_variant, am_pathogenicity, am_class)

        c_cv <- c_cv |>
            select(uniprot_id, protein_variant, transcript_id, cv_class)

        # Join databases by protein_variant
        res <-
            as.numeric(gsub(".*?([0-9]+).*", "\\1", c_am$protein_variant))

        c_am <- c_am |> mutate(aa_pos = res)

        c_combo <- left_join(c_am, c_cv,
                             by = c('uniprot_id', 'protein_variant')) |>
            relocate('transcript_id', .after = 'uniprot_id')

        # Check if protein has multiple transcripts
        unique_transcript_id <-
            c_combo |>
            filter(!is.na(transcript_id)) |>
            select(transcript_id) |>
            distinct() |>
            pull()

         if (length(unique_transcript_id) > 1L){
             stop("Multiple transcripts detected")
         }

        # Mutate the new column with protein transcript
        c_combo <- c_combo |>
            mutate(transcript_id = rep(unique_transcript_id))


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
            labs(title = paste0("UniProt ID: ", protein)) +
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
