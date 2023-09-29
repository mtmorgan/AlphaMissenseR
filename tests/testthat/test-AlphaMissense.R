test_that("ALPHAMISSENSE_RECORD is defined", {
    expect_identical(ALPHAMISSENSE_RECORD, "8360242")
})

test_that("'am_record_json()' returns", {
    json <- am_record_json(ALPHAMISSENSE_RECORD)
    expect_true(is_scalar_character(json))
})

test_that("'am_data_license()' reports appropriate license", {
    expect_identical(am_data_license(ALPHAMISSENSE_RECORD), "CC-BY-NC-SA-4.0")
})

test_that("'am_available()' works", {
    ## new cache
    fl <- tempfile(); dir.create(fl)
    bfc <- BiocFileCache::BiocFileCache(fl)
    db <- db_connect(ALPHAMISSENSE_RECORD, bfc)

    result <- am_available(ALPHAMISSENSE_RECORD, bfc)
    expect_true(inherits(result, "tbl_df"))
    expect_identical(NROW(result), 7L)
    expect_true(!any(result$cached))

    keys <- c(
        "aa_substitutions", "gene_hg19", "gene_hg38", "hg19",
        "hg38", "isoforms_aa_substitutions", "isoforms_hg38"
    )
    expect_true(setequal(keys, result$key))

    db_disconnect(db)
})

test_that("'am_data()' works", {
    ## hmm, not sure how to test without excessive download...
    expect_error(am_data("unknown_key"))

    fl <- tempfile(); dir.create(fl)
    bfc <- BiocFileCache::BiocFileCache(fl)
    db <- db_connect(ALPHAMISSENSE_RECORD, bfc)
    ## should create empty database before failing
    expect_error(am_data("unknown_key", bfc))
    expect_identical(BiocFileCache::bfccount(bfc), 1L)
    expect_identical(db_tables(db), character(0))

    db_disconnect(db)
})
