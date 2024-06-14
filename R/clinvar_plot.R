#' Filter the AlphaMissense table with uniprotID
#'
#' @noRd
#'
#' @importFrom dplyr filter as_tibble
#'
clinvar_filter_am_table <-
    function(am_table, uID)
{
    ## Check if am_table is missing
    if (missing(am_table)) {
        spdl::info(paste0(
            "'alphamissense_table' not provided, using default ",
            "'am_data(\"aa_substitution\")' table accessed through ",
            "the AlphaMissenseR package"
        ))

        am_table <- am_data("aa_substitutions")
    }

    ## Take alphamissense_table and filter for the uniprotId
    alphamissense_table <-
        am_table |>
        filter(.data$uniprot_id == uID) |>
        as_tibble()

    ## Check if table is empty after filtering
    ## This will work for a tibble or a data.frame
    if (!NROW(alphamissense_table)) {
        stop(
            "no AlphaMissense information found for the protein ",
            "accession '", uID, "'; check that the UniProt ID is correct"
        )
    }

    alphamissense_table
}


#' Filter the clinvar table with uniprot ID
#'
#' @noRd
#'
#' @importFrom dplyr filter
#'
clinvar_filter_cv_table <-
    function(cv_table, uID)
{
    if (missing(cv_table)) {
        spdl::info(paste0(
            "'clinvar_table' not provided, using default ",
            "ClinVar dataset in AlphaMissenseR package"
        ))

        clinvar <- clinvar_data()

        # Separate UniProt ID and Protein Variant
        cv_table <-
            clinvar |>
            as_tibble() |>
            tidyr::separate(
                .data$protein_variant,
                into = c("uniprot_id", "protein_variant"),
                sep = ":"
            ) |>
            mutate(
                cv_class = as.factor(.data$label)
            )
    }

    ## Take clinvar_table and filter for the uniprotId
    clinvar_table <-
        cv_table |>
        filter(.data$uniprot_id == uID)


    ## Check if the table is empty after filtering
    if (!NROW(clinvar_table)) {
        stop(
            "no ClinVar information found for the protein ",
            "accession '", uID, "'; check that the UniProt ID is correct"
        )
    }

    clinvar_table
}

#' Prepare data for the function clinvar_plot
#'
#' @noRd
#'
#' @importFrom dplyr left_join mutate case_when mutate_at group_by
#'     ungroup arrange
#'
clinvar_prepare_data_for_plot <-
    function(am_table, cv_table)
{
    ## grab amino acid positions
    am_table <-
        am_table |>
        mutate(
            aa_pos = as.integer(
                gsub(".*?([0-9]+).*", "\\1", .data$protein_variant)
            ))

    ## join datasets
    combined_data <-
        left_join(
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
                    "AM ambiguous")
        ) |>
        mutate_at(vars(.data$code_color), factor) |>
        arrange(.data$code_color)

    ## Grab the thresholds for AM pathogenicity to plot
    combined_data <-
        combined_data |>
        group_by(.data$am_class) |>
        mutate(
            max = max(.data$am_pathogenicity, na.rm=TRUE),
            min = min(.data$am_pathogenicity, na.rm=TRUE)
        ) |>
        ungroup()

    combined_data
}

#' Create a ClinVar plotting function using ggplot
#'
#' @noRd
#'
#' @importFrom ggplot2 ggplot geom_point aes scale_colour_manual element_text
#'     scale_fill_manual scale_shape_manual scale_size_manual element_blank
#'     scale_discrete_manual geom_hline labs xlab ylab theme_classic
#'     theme ggplot_build
#'
clinvar_create_plot <-
    function(combined_table, uId)
{
    ## Create named vectors for all scale layers
    colScale <-
        scale_colour_manual(
        name = "code_color",
        values = c(
            "AM ambiguous" = "gray",
            "AM benign" = "#89d5f5",
            "AM pathogenic" = "#f56c6c",
            "CV benign" = "black",
            "CV pathogenic" = "black"
            ))

    fillScale <-
        scale_fill_manual(
        name = "code_color",
        values = c(
            "AM ambiguous" = "gray",
            "AM benign" = "#89d5f5",
            "AM pathogenic" = "#f56c6c",
            "CV benign" = "#007cb0",
            "CV pathogenic" = "#c70606"
            ))

    shapeScale <-
        scale_shape_manual(
        name = "code_color",
        values = c(
            "AM ambiguous" = 19,
            "AM benign" = 19,
            "AM pathogenic" = 19,
            "CV benign" = 21,
            "CV pathogenic" = 21
            ))

    sizeScale <-
        scale_size_manual(
        name = "code_color",
        values = c(
            "AM ambiguous" = 2,
            "AM benign" = 2,
            "AM pathogenic" = 2,
            "CV benign" = 4,
            "CV pathogenic" = 4
            ))

    strokeScale <-
        scale_discrete_manual(
        name = "code_color",
        aesthetics = "stroke",
        values = c(
            "AM ambiguous" = 0,
            "AM benign" = 0,
            "AM pathogenic" = 0,
            "CV benign" = 1.5,
            "CV pathogenic" = 1.5
            ))

    cv_plot <-
        combined_table |>
        ggplot(aes(.data$aa_pos, .data$am_pathogenicity)) +
        geom_point(aes(
            shape = .data$code_color,
            color = .data$code_color,
            size = .data$code_color,
            fill = .data$code_color,
            stroke = .data$code_color
        )) +
        strokeScale +
        shapeScale +
        sizeScale +
        colScale +
        fillScale +
        geom_hline(
            yintercept =
                combined_table |>
                filter(.data$am_class == "pathogenic") |>
                pull(min) |>
                unique(),
            linetype = 2,
            color = "#c70606"
        ) +
        geom_hline(
            yintercept =
                combined_table |>
                filter(.data$am_class == "benign") |>
                pull(max) |>
                unique(),
            linetype = 2,
            color = "#007cb0"
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

#' @rdname clinvar_plot
#'
#' @title Integrate ClinVar Labels with AlphaMissense Pathogenicity Scores
#'
#' @description `clinvar_plot()` integrates ClinVar classifications
#'    with AlphaMissense predicted scores derived from
#'    `am_data("aa_substitutions")` and returns a ggplot object for
#'    visualization.
#'
#' @param uniprotId `character()` a valid UniProt accession identifier.
#'
#' @param alphamissense_table a table containing AlphaMissense
#'    predictions for protein variants. By default, the table is
#'    derived from `am_data("aa_substitution")`. Alternatively, a
#'    user-defined [`tibble::tbl_df`] or [`data.frame`]
#'    can be supplied.
#'
#' @param clinvar_table a table containing ClinVar information. By
#'    default, the table is derived from the supplemental data of the
#'    AlphaMissense paper. Alternatively, a user-defined
#'    [`tibble::tbl_df`] or [`data.frame`] can be
#'    supplied.
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
#' @return `clinvar_plot()` returns a `ggplot` object which overlays
#'    ClinVar classifications onto AlphaMissense predicted
#'    scores. Blue, gray, and red colors represent pathogenicity
#'    classifications for "likely benign", "ambiguous", or
#'    "likely pathogenic", respectively. Large, bolded points are
#'    ClinVar variants colored according to their clinical
#'    classification, while smaller points in the background are
#'    AlphaMissense predictions.
#'
#' @examples
#'
#' alphamissense_table <- am_data("aa_substitutions")
#'
#' clinvar_plot(uniprotId = "P37023",
#'    alphamissense_table = alphamissense_table)
#'
#' @references Cheng et al.,
#' Accurate proteome-wide missense variant effect prediction with AlphaMissense.
#' \emph{Science} 381, eadg7492. DOI:10.1126/science.adg7492.
#'
#' @importFrom BiocBaseUtils isCharacter
#'
#' @importFrom utils data
#'
#' @export
clinvar_plot <-
    function(uniprotId, alphamissense_table, clinvar_table)
{
    ## Validate uniprotId to start filtering alphamissense and clinvar tables
    stopifnot(isCharacter(uniprotId))

    ## Filter AM and CV tables with uniProtID
    alphamissense_table <-
        clinvar_filter_am_table(
            am_table = alphamissense_table,
            uID = uniprotId
        )

    clinvar_table <-
        clinvar_filter_cv_table(
            cv_table = clinvar_table,
            uID = uniprotId
        )

    ## Validate extracted tables
    stopifnot(
        is.data.frame(alphamissense_table) | is.tbl(alphamissense_table),
        all(c("uniprot_id", "protein_variant", "am_class", "am_pathogenicity")
            %in% colnames(alphamissense_table)),
        is.data.frame(clinvar_table) | is.tbl(clinvar_table),
        all(c("uniprot_id", "protein_variant", "cv_class")
            %in% colnames(clinvar_table))
    )

    ## Process tables for plotting
    combined_table <-
        clinvar_prepare_data_for_plot(
            am_table = alphamissense_table,
            cv_table = clinvar_table
        )

    ## Plot sequence window
    clinvar_create_plot(combined_table = combined_table, uId = uniprotId)
}

#' @rdname clinvar_plot
#'
#' @description `clinvar_data()` loads in the raw ClinVar information from
#'    the supplemental table of the AlphaMissense
#'    [\[2023\]](https://www.science.org/doi/10.1126/science.adg7492) paper.
#'
#' @return
#'
#' `clinvar_data()` returns a tbl with 82872 rows and 5 variables:
#'
#' - `variant_id`: ClinVar variant identifer.
#' - `transcript_id`: Ensembl transcript identifier.
#' - `protein_variant`: UniProt accession:protein variant identifier.
#' - `AlphaMissense`: AlphaMissense pathogenicity score.
#' - `label`: Binary ClinVar class. 0 for benign or 1 for pathogenic.
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
