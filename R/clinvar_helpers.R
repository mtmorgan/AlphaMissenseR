# plot_clinvar.R helper functions

.check_am_table <-
    function(am_table, uID){
        if (missing(am_table)) {
        message("'alphamissense_table' not provided, using default ",
                "'am_data(\"aa_substitution\")' table accessed through ",
                "the AlphaMissenseR package.")

        alphamissense_table <- am_data("aa_substitutions") |>
            filter(.data$uniprot_id == uID) |>
            dplyr::collect()

        } else {

        # Take user-defined alphamissense_table and filter for the uniprotId
        alphamissense_table <- am_table |>
            filter(.data$uniprot_id == uID)
        }

        return(alphamissense_table)
}


.check_cv_table <-
    function(cv_table, uID) {

        if (missing(cv_table)) {
            message("'clinvar_table' not provided, using default ",
                "ClinVar dataset in AlphaMissenseR package")

        data_env <- new.env(parent = emptyenv())
        data("clinvar_data", envir = data_env, package = "AlphaMissenseR")
        clinvar_data <- data_env[["clinvar_data"]]

        clinvar_table <-  clinvar_data |>
            filter(.data$uniprot_id == uID)

    } else {
        # Take user-defined clinvar_table and filter for the uniprotId
        clinvar_table <- cv_table |>
            filter(.data$uniprot_id == uID)
    }

    return(clinvar_table)
}


.process_data_for_plot_clinvar <-
    function(am_table = alphamissense_table, cv_table = clinvar_table)
{
        # grab amino acid positions
        am_table <- am_table |>
            mutate(
                aa_pos = as.numeric(
                    gsub(".*?([0-9]+).*", "\\1", .data$protein_variant)
                )
            )

        # join datasets
        c_combo <- left_join(
            am_table,
            cv_table,
            by = c('uniprot_id', 'protein_variant')
        )

        # add color code matching AM and CV labels
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

    return(c_combo)
}

.make_cv_plot <-
    function(combo_table, cutoff, uId)
{
    cv_plot <- combo_table |>
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
                        yintercept = cutoff$path_cutoff, linetype = 2,
                        color = "#c70606"
                    ) +
                    geom_hline(
                        yintercept = cutoff$benign_cutoff,
                        linetype = 2, color = "#007cb0"
                    ) +
                    labs(title = paste0("UniProt ID: ", uId)) +
                    xlab("amino acid position") +
                    ylab("AlphaMissense score") +
                    theme_classic()

                cv_plot +
                    theme(
                        axis.text.x = element_text(size = 16),
                        axis.text.y = element_text(size = 16),
                        axis.title.y = element_text(size = 16),
                        axis.title.x = element_text(size = 16),
                        legend.title = element_blank(),
                        legend.text = element_text(size = 11)
                    )
            }
