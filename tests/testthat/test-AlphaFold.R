test_that("AlphaFold endpoints exist", {
    skip_if_offline()

    ## af_predictions()
    url <- whisker.render(ALPHAFOLD_PREDICTION, c(qualifier = "P35557"))
    ## 'HEAD()' returns 404
    expect_identical(httr::status_code(httr::GET(url)), 200L)

    ## af_prediction_view() -- PDB file
    url <- "https://alphafold.ebi.ac.uk/files/AF-P35557-F1-model_v4.pdb"
    expect_identical(httr::status_code(httr::HEAD(url)), 200L)
})
