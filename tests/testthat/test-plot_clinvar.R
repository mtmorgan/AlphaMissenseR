test_that("prepare_data_for_plot_clinvar() works", {
    ## All the data munging happens within this function
    am_data <- data.frame(uniprot_id = c("A", "A", "A", "A", "A", "B"),
                          protein_variant = c("M1A", "M2C",
                                              "R3L", "K4L", "G5L", "F50G"),
                          am_pathogenicity = c(0.25, 0.89, 0.5, 0.2, 0.40, 0.60),
                          am_class = c("benign", "pathogenic", "ambiguous",
                                       "benign", "ambiguous", "pathogenic"))

    ## TODO: Add demo data for cv_class = 1
    cv_data <- data.frame(uniprot_id = c("A", "A", "A", "A"),
                      protein_variant = c("M1A", "M2C", "G5L", "P10R"),
                      cv_class = c(0, 0, 0, 0))

    res <- prepare_data_for_plot_clinvar(
        am_table = am_data,
        cv_table = cv_data
    )

    col_names <- c("uniprot_id", "protein_variant", "am_pathogenicity",
                   "am_class", "aa_pos", "cv_class", "code_color",
                   "max", "min")

    ## All colnames are correctly present
    expect_true(all(names(res) %in% col_names))

    ## There are only 3 am_class types - ambiguous, benign, pathogenic
    expect_equal(res |> distinct(am_class) |> nrow(), 3L)

    ## TODO: check code_color in the res and make they are valid

    ## TODO: check max and min for pathogenecity

    ## TODO: check that returned data frame is "grouped"
})


## TODO: Fill these up
test_that("filter_am_table() works", {

    ## Write test case when am_table is missing

    ## Write test case when bad am_table is given i.e
    ## 1. uniprotID with no hits

    ## Write test case when am_table is user provided
    ## ex: filter to a subset of the aa_substitutions table

})


## TODO: Fill these up
test_that("filter_cv_table() works", {

    ## Write test case when cv_table is missing

    ## Write test case when bad cv_table is given i.e
    ## 1. uniprot ID with no hits

    ## Write test case when cv_table is user provided
    ## ex: filter to a subset of the clinvar data table

})



test_that("'plot_clinvar()' works", {
    ## Create example dataframes for alphamissense and clinvar data
    df1 <- data.frame(uniprot_id = c("A", "A", "A", "A", "A", "B"),
                      protein_variant = c("M1A", "M2C",
                                          "R3L", "K4L", "G5L", "F50G"),
                      am_pathogenicity = c(0.25, 0.89, 0.5, 0.2, 0.40, 0.60),
                      am_class = c("benign", "pathogenic", "ambiguous",
                                   "benign", "ambiguous", "pathogenic"))

    df2 <- data.frame(uniprot_id = c("A", "A", "A", "A"),
                      protein_variant = c("M1A", "M2C", "G5L", "P10R"),
                      cv_class = c(0, 0, 0, 0))

    # Call the function with example input
    plot <- plot_clinvar(uniprotId = "A",
                         alphamissense_table = df1,
                         clinvar_table = df2)

    # check plot axes
    expect_identical(plot$labels$y, "AlphaMissense score")
    expect_identical(plot$labels$x, "amino acid position")

    # check plot layers
    expect_in(class(plot$layers[[1]]$geom),
              c("GeomPoint", "Geom", "ggproto", "gg"))
    expect_in(class(plot$layers[[2]]$geom),
              c("GeomHline", "Geom", "ggproto", "gg"))
    expect_in(class(plot$layers[[3]]$geom),
              c("GeomHline", "Geom", "ggproto", "gg"))

    # check plot scales
    expect_in(plot$scales$scales[[1]]$aesthetics, "stroke")
    expect_in(plot$scales$scales[[2]]$aesthetics, "shape")
    expect_in(plot$scales$scales[[3]]$aesthetics, "size")
    expect_in(plot$scales$scales[[4]]$aesthetics, "colour")
    expect_in(plot$scales$scales[[5]]$aesthetics, "fill")

    # check if the result is a ggplot object
    expect_true("ggplot" %in% class(plot))

    ## TODO:
    ## 1. Write test case where numeric uniprot ID is given
    expect_error(plot_clinvar(uniprotId = "A"))

    ## TODO
    ## 2. Write a test case where AM table with bad colnames is given

    ## TODO
    ## 3. Write test case where CV table with bad colnames is given

})
