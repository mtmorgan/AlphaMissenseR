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

})
