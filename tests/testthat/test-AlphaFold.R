test_that("AlphaFold endpoints exist", {
    ## skip_if_offline() also checks for environment variable 'NOT_CRAN'
    skip_if_not_installed("curl")
    skip_if_not(!is.null(curl::nslookup("captive.apple.com")))

    ## af_predictions()
    url <- whisker.render(ALPHAFOLD_PREDICTION, c(qualifier = "P35557"))
    ## 'HEAD()' returns 404
    expect_identical(httr::status_code(httr::GET(url)), 200L)

    ## af_prediction_view() -- PDB file
    url <- "https://alphafold.ebi.ac.uk/files/AF-P35557-F1-model_v4.pdb"
    expect_identical(httr::status_code(httr::HEAD(url)), 200L)
})

test_that("af_predictions() works", {
    ## skip_if_offline() also checks for environment variable 'NOT_CRAN'
    skip_if_not_installed("curl")
    skip_if_not(!is.null(curl::nslookup("captive.apple.com")))

    expect_error(
        af_predictions(character()),
        "length(uniprot_ids) >= 1L is not TRUE", fixed = TRUE
    )

    expect_output(
        tbl <- af_predictions(c("P35557", "xyz")),
        "1 of 2 uniprot accessions not found\n  'xyz'"
    )
    expect_true(NROW(tbl) == 1L)
    expect_true(NCOL(tbl) >= 23L)

    ## as of 22 May, 2024
    colnames <- c(
        "entryId", "gene", "uniprotAccession", "uniprotId",
        "uniprotDescription", "taxId", "organismScientificName",
        "uniprotStart", "uniprotEnd", "uniprotSequence",
        "modelCreatedDate", "latestVersion", "allVersions",
        "isReviewed", "isReferenceProteome", "cifUrl", "bcifUrl",
        "pdbUrl", "paeImageUrl", "paeDocUrl", "amAnnotationsUrl",
        "amAnnotationsHg19Url", "amAnnotationsHg38Url"
    )
    expect_true(all(colnames %in% names(tbl)))
})

test_that("af_colorfunc_by_position() works", {
    palette <- as.character(1:10) # proxy  palette, easier to test
    ## 20 positions, colored in pairs
    colors <- af_colorfunc_by_position_colors(
        pos = 1:20, value = 20:1 / 20, pos_max = 20L,
        palette = palette, palette_min = NULL, palette_max = NULL
    )
    expect_identical(colors, as.character(rep(10:1, each = 2)))

    ## 12 'pos'itions, no color information for position 1 or 12
    colors <- af_colorfunc_by_position_colors(
        pos = 1L + 1:10, value = 10:1 / 10, pos_max = 12L,
        palette = palette, palette_min = NULL, palette_max = NULL
    )
    expect_identical(colors, c("gray", 10:1, "gray"))

    ## 20 positions, compressed palette
    colors <- af_colorfunc_by_position_colors(
        pos = 1:20, value = 20:1 / 100, pos_max = 20L,
        palette = palette, palette_min = 0L, palette_max = 1L
    )
    expect_identical(colors, as.character(rep(2:1, each = 10)))
})

test_that("af_colorfunc_by_position() works", {
    ## see test 'af_colorfunc_by_position_colors()' for test of color
    ## mapping
    df <- tibble(
        pos = 1 + 1:10, # no color information for position 1
        value = 10:1 / 10
    )

    ## hard to test validity of returnd javascript function
    expect_identical(nchar(af_colorfunc_by_position(df, "pos", "value")), 205L)
})

test_that("af_prediction_view() returns an appropriate object", {
    skip_if_not_installed("curl")
    skip_if_not(!is.null(curl::nslookup("captive.apple.com")))

    skip_if_not_installed("r3dmol") # Suggests: packages
    skip_if_not_installed("bio3d")

    P35557 <-
        am_data("hg38") |>
        dplyr::filter(uniprot_id == "P35557")
    view <- af_prediction_view(P35557)

    ## hard to test details of returned object
    expect_s3_class(view, "r3dmol")

    P35555 <-
        am_data("hg38") |>
        dplyr::filter(uniprot_id == "P35555")
    ## expect_error() fails to detect error, for some reason
    err <- tryCatch(af_prediction_view(P35555), error = identity)
    expect_s3_class(err, "simpleError")
    expect_identical(
        conditionMessage(err), 
        "'af_prediction_view()' could not find UniProt accession 'P35555'"
    )
})
