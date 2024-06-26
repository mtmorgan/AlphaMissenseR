test_that("clinvar_filter_am_table() works", {
    ## Test case when am_table is missing
    spdl::set_level("warn")
    on.exit(spdl::set_level("info"))
    res <- clinvar_filter_am_table(uID = "P37023")

    expect_equal(res |> NROW(), 9557L)

    ## Test case when am_table is provided by user
    am_data <- data.frame(
        uniprot_id = c("A", "A", "A", "A", "A", "B"),
        protein_variant = c("M1A", "M2C", "R3L", "K4L", "G5L", "F50G"),
        am_pathogenicity = c(0.2, 0.8, 0.5, 0.2, 0.4, 0.6),
        am_class = c(
            "benign", "pathogenic", "ambiguous", "benign",
            "ambiguous", "pathogenic"
        )
    )
    object <- clinvar_filter_am_table(uID = "A", am_data)

    expected <- tibble(
        uniprot_id = c("A", "A", "A", "A", "A"),
        protein_variant = c("M1A", "M2C", "R3L", "K4L", "G5L"),
        am_pathogenicity = c(0.2, 0.8, 0.5, 0.2, 0.4),
        am_class = c(
            "benign", "pathogenic", "ambiguous", "benign",
            "ambiguous"
        )
    )

    expect_identical(object, expected)

    ## Test case when bad am_table is given i.e uniprotID with no hits
    expect_error(
        clinvar_filter_am_table(uID = "C", am_table = am_data),
        paste(
            "no AlphaMissense information found for the protein accession",
            "'C'; check that the UniProt ID is correct"
        )
    )
})

test_that("clinvar_filter_cv_table() works", {
    ## Test case when cv_table is missing
    spdl::set_level("warn")
    on.exit(spdl::set_level("info"))
    res <- clinvar_filter_cv_table(uID = "P37023")

    expect_equal(res |> NROW(), 113L)

    ## Test case when cv_table is provided by user
    cv_data <- data.frame(
        uniprot_id = c("A", "A", "A", "A", "B"),
        protein_variant = c("M1A", "M2C", "G5L", "P10R", "F50G"),
        cv_class = c("pathogenic", rep("benign", 3L), "pathogenic")
    )
    object <- clinvar_filter_cv_table(uID = "A", cv_table = cv_data)

    expected <- tibble(
        uniprot_id = c("A", "A", "A", "A"),
        protein_variant = c("M1A", "M2C", "G5L", "P10R"),
        cv_class = c("pathogenic", rep("benign", 3L))
    )

    expect_identical(object, expected)

    ## Test case when bad cv_table is given i.e uniprotID with no hits
    expect_error(
        clinvar_filter_cv_table(uID = "C", cv_table = cv_data),
        paste(
            "no ClinVar information found for the protein",
            "accession 'C'; check that the UniProt ID is correct"
        )
    )
})

test_that("clinvar_prepare_data_for_plot() works", {
    ## All the data munging happens within this function
    am_data <- data.frame(
        uniprot_id = c("A", "A", "A", "A", "A", "B"),
        protein_variant = c("M1A", "M2C", "R3L", "K4L", "G5L", "F50G"),
        am_pathogenicity = c(0.2, 0.8, 0.5, 0.2, 0.4, 0.6),
        am_class = c(
            "benign", "pathogenic", "ambiguous", "benign", "ambiguous",
            "pathogenic"
        )
    )

    cv_data <- data.frame(
        uniprot_id = c("A", "A", "A", "A"),
        protein_variant = c("M1A", "M2C", "G5L", "P10R"),
        cv_class = c("pathogenic", rep("benign", 3L))
    )

    object <- suppressWarnings(clinvar_prepare_data_for_plot(
        am_table = am_data,
        cv_table = cv_data
    ))

    expected <- tibble(
        uniprot_id = c("A", "A", "B", "A", "A", "A"),
        protein_variant = c("R3L", "K4L", "F50G", "M2C", "G5L", "M1A"),
        am_pathogenicity = c(0.5, 0.2, 0.6, 0.8, 0.4, 0.2),
        am_class = c("ambiguous", "benign", "pathogenic", "pathogenic",
                     "ambiguous", "benign"),
        aa_pos = c(3L, 4L, 50L, 2L, 5L, 1L),
        cv_class = c(NA, NA, NA, "benign", "benign", "pathogenic"),
        code_color = as.factor(c(
            "AM ambiguous", "AM benign", "AM pathogenic", "CV benign",
            "CV benign", "CV pathogenic"
        )),
        max = c(0.5, 0.2, 0.8, 0.8, 0.5, 0.2),
        min = c(0.4, 0.2, 0.6, 0.6, 0.4, 0.2)
    )

    expect_identical(object, expected)

    ## All colnames are correctly present
    col_names <- c(
        "uniprot_id", "protein_variant", "am_pathogenicity",
        "am_class", "aa_pos", "cv_class", "code_color", "max", "min"
    )

    expect_true(all(names(object) %in% col_names))

    ## There are only 3 am_class types - ambiguous, benign, pathogenic
    expect_equal(object |> distinct(am_class) |> NROW(), 3L)

    ## There are only 5 code_color types - check assigned correctly
    color_names <- c(
        "AM ambiguous", "AM benign", "AM pathogenic", "CV benign",
        "CV benign", "CV pathogenic"
    )

    expect_identical(
        object |>
            select(code_color) |>
            pull() |>
            as.vector(),
        color_names
    )

    ## Check max and min for pathogenicity were calculated correctly
    path_scores <- c(0.6, 0.8, 0.2, 0.2)

    expect_identical(
        object |>
            filter(am_class == "pathogenic") |>
            pull(min) |>
            unique(),
        path_scores[1]
    )
    expect_identical(
        object |>
            filter(am_class == "pathogenic") |>
            pull(max) |>
            unique(),
        path_scores[2]
    )
    expect_identical(
        object |>
            filter(am_class == "benign") |>
            pull(min) |>
            unique(),
        path_scores[3]
    )
    expect_identical(
        object |>
            filter(am_class == "benign") |>
            pull(max) |>
            unique(),
        path_scores[4]
    )
})

test_that("'clinvar_plot()' works", {
    ## Create example dataframes for alphamissense and clinvar data
    am_data <- data.frame(
        uniprot_id = c("A", "A", "A", "A", "A", "B"),
        protein_variant = c("M1A", "M2C", "R3L", "K4L", "G5L", "F50G"),
        am_pathogenicity = c(0.2, 0.8, 0.5, 0.2, 0.4, 0.6),
        am_class = c(
            "benign", "pathogenic", "ambiguous", "benign", "ambiguous",
            "pathogenic"
        )
    )

    cv_data <- data.frame(
        uniprot_id = c("A", "A", "A", "A"),
        protein_variant = c("M1A", "M2C", "G5L", "P10R"),
        cv_class = c("pathogenic", rep("benign", 3))
    )

    ## Call the function with example input
    plot <- suppressWarnings(clinvar_plot(
        uniprotId = "A",
        alphamissense_table = am_data,
        clinvar_table = cv_data
    ))

    ## Check aesthetics values are assigned correctly to groups
    stroke_vec <- c(0.0, 0.0, 1.5, 1.5, 1.5)
    expect_identical(ggplot_build(plot)$data[[1]]$stroke, stroke_vec)

    shape_vec <- c(19, 19, 21, 21, 21)
    expect_identical(ggplot_build(plot)$data[[1]]$shape, shape_vec)

    size_vec <- c(2, 2, 4, 4, 4)
    expect_identical(ggplot_build(plot)$data[[1]]$size, size_vec)

    color_vec <- c("gray", "#89d5f5", "black", "black", "black")
    expect_identical(ggplot_build(plot)$data[[1]]$colour, color_vec)

    fill_vec <- c("gray", "#89d5f5", "#007cb0", "#007cb0", "#c70606")
    expect_identical(ggplot_build(plot)$data[[1]]$fill, fill_vec)

    ## check plot axes
    expect_identical(plot$labels$y, "AlphaMissense score")
    expect_identical(plot$labels$x, "amino acid position")

    ## check plot layers
    expect_in(
        class(plot$layers[[1]]$geom),
        c("GeomPoint", "Geom", "ggproto", "gg")
    )
    expect_in(
        class(plot$layers[[2]]$geom),
        c("GeomHline", "Geom", "ggproto", "gg")
    )
    expect_in(
        class(plot$layers[[3]]$geom),
        c("GeomHline", "Geom", "ggproto", "gg")
    )

    ## check plot scales exists in aesthetics
    expect_in(plot$scales$scales[[1]]$aesthetics, "stroke")
    expect_in(plot$scales$scales[[2]]$aesthetics, "shape")
    expect_in(plot$scales$scales[[3]]$aesthetics, "size")
    expect_in(plot$scales$scales[[4]]$aesthetics, "colour")
    expect_in(plot$scales$scales[[5]]$aesthetics, "fill")

    ## check if the result is a ggplot object
    expect_true("ggplot" %in% class(plot))
})
