#' @importFrom BiocBaseUtils isCharacter
#' @importFrom ggplot2 ggplot
#'
#' @export
plot_clinvar <-
    function(protein, clinvar)
    {
        # validity checks
        stopifnot(isCharacter(protein))

        protein <- as.character(protein)

        c_cv <- clinvar |>
            filter(accession == protein)

        # Check if protein found in ClinVar
        if (nrow(c_cv) < 1)
            stop(
                "No ClinVar information found for the protein accession. ",
                "Check that the UniProt ID is correct."
            )

        # Grab AM proteins that match CV
        c_AM <- db_connect() |>
            tbl("aa_substitutions") |>
            filter(uniprot_id == protein) |>
            dplyr::as_tibble()

        # Check if protein found in AM
        if (nrow(c_AM) < 1) {
            stop(
                "No AlphaMissense information found for the protein accession.",
                "Check that the UniProt ID is correct."
            )
        }

        # Join databases by protein variant
        res <-
            as.numeric(gsub(".*?([0-9]+).*", "\\1", c_AM$protein_variant))
        c_AM <- c_AM |> mutate(aa_pos = res)

        c_combo <- left_join(c_AM, c_cv, by = "protein_variant") |>
            select(-c('accession', 'AlphaMissense')) |>
            rename(cv_variant_id = variant_id, cv_label = label) |>
            relocate('transcript_id', .after = 'uniprot_id')

        # Check if protein has multiple transcripts
        if (length(unique(c_combo$transcript_id[!is.na(c_combo$transcript_id)])) > 1)
            stop(c("Multiple transcripts detected."))

        # change to tidyverse way
        # make a new variable
        # unique_transcript_id <-
        #     c_combo |>
        #     filter(!is.na(transcript_id)) |>
        #     distinct()
        # if (length(unique_transcript_id) > 1L)
        #     stop("Multiple transcripts detected")

        # change to tidyverse way to mutate the new column
        c_combo$transcript_id <-
            rep(unique(c_combo$transcript_id[!is.na(c_combo$transcript_id)]))

        # c_combo <- c_combo |>
        #     mutate(transcript_id = rep(unique_transcript_id))


        ##### PLOTTING #####

        # Grab the thresholds for pathogenicity/benign to plot (dashed lines)
        am_cutoff <- c_combo |>
            filter(am_class == "ambiguous") |>
            select(am_pathogenicity) |>
            summarise(
                path_cutoff = min(am_pathogenicity),
                benign_cutoff = max(am_pathogenicity)
            )


        # Add color code matching AM and CV labels
        c_combo <- c_combo |>
            mutate(
                code_color = case_when(
                    !is.na(cv_label) & cv_label == "0" ~ "CV_benign",
                    !is.na(cv_label) &
                        cv_label == "1" ~ "CV_pathogenic",
                    is.na(cv_label) &
                        am_class == "pathogenic" ~ "AM_pathogenic",
                    is.na(cv_label) &
                        am_class == "benign" ~ "AM_benign",
                    is.na(cv_label) &
                        am_class == "ambiguous" ~ "AM_ambiguous"
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
