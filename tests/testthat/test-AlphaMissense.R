test_that("ALPHAMISSENSE_RECORD is defined", {
    expect_identical(ALPHAMISSENSE_RECORD, "8360242")
})

test_that("'am_record_json()' returns", {
    json <- am_record_json(ALPHAMISSENSE_RECORD)
    expect_true(is_scalar_character(json))
})

test_that("'am_data_license()' reports appropriate license", {
    output <- capture.output({
        result <- am_data_license(ALPHAMISSENSE_RECORD)
    })
    expect_identical(result, "CC-BY-NC-SA-4.0")
    expect_true(nzchar(output))
})

test_that("'am_available()' works", {
    ## new cache
    fl <- tempfile(); dir.create(fl)
    bfc <- BiocFileCache::BiocFileCache(fl)

    result <- am_available(ALPHAMISSENSE_RECORD, bfc)
    expect_true(inherits(result, "tbl_df"))
    expect_identical(NROW(result), 7L)
    expect_true(!any(result$cached))

    keys <- c(
        "aa_substitutions", "gene_hg19", "gene_hg38", "hg19",
        "hg38", "isoforms_aa_substitutions", "isoforms_hg38"
    )
    expect_true(setequal(keys, result$key))
})

test_that("'am_data_import_csv()' works", {
    record <- ALPHAMISSENSE_RECORD
    fl <- tempfile(); dir.create(fl)
    bfc <- BiocFileCache::BiocFileCache(fl)

    tsv_file <- tempfile()
    writeLines(c(
        "# Comment",
        "#",
        "# Comment",
        "#CHROM\tPOS",
        "chr1\t1",
        "chr1\t2"
    ), tsv_file)

    spdl::set_level("warn")
    output <- capture.output({
        tbl <- am_data_import_csv(record, bfc, "tsv_file", tsv_file)
    }, type = "message")
    spdl::set_level("info")
    expect_true(nzchar(output))
    expect_true(NROW(tbl |> collect()) == 2L)
    expect_identical(colnames(tbl), c("CHROM", "POS"))

    db_disconnect(db_connect(record, bfc))
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
