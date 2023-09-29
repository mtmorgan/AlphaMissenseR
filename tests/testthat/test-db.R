test_that("'db_connect()' works", {
    fl <- tempfile(); dir.create(fl)
    bfc <- BiocFileCache::BiocFileCache(fl)

    ## read-only connection
    db <- db_connect(ALPHAMISSENSE_RECORD, bfc)
    expect_true(DBI::dbIsValid(db))
    expect_identical(db, db_connect(ALPHAMISSENSE_RECORD, bfc))

    ## read-write connection
    db_rw <- db_connect(ALPHAMISSENSE_RECORD, bfc, read_only = FALSE)
    expect_true(DBI::dbIsValid(db))
    expect_identical(
        db_rw, db_connect(ALPHAMISSENSE_RECORD, bfc, read_only = FALSE)
    )

    ## different connections
    expect_true(!identical(db, db_rw))
    db_unmanaged <- db_connect(ALPHAMISSENSE_RECORD, bfc, managed = FALSE)
    expect_true(DBI::dbIsValid(db_unmanaged))
    expect_true(!identical(db, db_unmanaged))
    db_disconnect(db_unmanaged)

    db_disconnect(db_rw)
    db_disconnect(db)
})

test_that("'db_tables() works", {
    fl <- tempfile(); dir.create(fl)
    bfc <- BiocFileCache::BiocFileCache(fl)

    db <- db_connect(ALPHAMISSENSE_RECORD, bfc)
    expect_identical(db_tables(db), character(0))
    db_disconnect(db)

    ## new table seen by db_table(); requires rw connection
    db_rw <- db_connect(ALPHAMISSENSE_RECORD, bfc, read_only = FALSE)
    expect_identical(db_tables(db_rw), character(0))
    dbWriteTable(db_rw, "mtcars", mtcars)
    expect_identical(db_tables(db_rw), "mtcars")
    db_disconnect(db_rw)

    ## table also seen by (new) read-only connection
    db <- db_connect(ALPHAMISSENSE_RECORD, bfc)
    expect_identical(db_tables(db), "mtcars")
    db_disconnect(db)
})

test_that("'db_temporary_table()' works", {
    fl <- tempfile(); dir.create(fl)
    bfc <- BiocFileCache::BiocFileCache(fl)

    ## can't write table to read-only database
    db <- db_connect(ALPHAMISSENSE_RECORD, bfc)
    expect_error(dbWriteTable(db, "mtcars", mtcars))

    ## can write temporary table
    tbl <- db_temporary_table(db, mtcars, "mtcars")
    expect_identical(NROW(mtcars), count(tbl) |> pull(n) |> as.integer())

    db_disconnect(db)
})

test_that("'db_range_join()' works", {
    fl <- tempfile(); dir.create(fl)
    bfc <- BiocFileCache::BiocFileCache(fl)

    db <- db_connect(ALPHAMISSENSE_RECORD, bfc)

    key <- tibble(
        `#CHROM` = paste0("chr", rep(1:2, c(3,5))),
        POS = 10 * c(1:3, 1:5),
        key_etc = letters[1:8]
    )
    key_dbtbl <- db_temporary_table(db, key, "key")

    ## empty join
    join <- tibble(
        `#CHROM` = character(),
        start = integer(),
        end = integer(),
        join_etc = character()
    )
    join_dbtbl <- db_temporary_table(db, join, "join_tbl0")
    result <- db_range_join(db, "key", "join_tbl0", "res0") |> collect()
    expect_identical(NROW(result), 0L)
    expect_true(setequal(
        colnames(result),
        unique(c(colnames(key), colnames(join)))
    ))

    ## single join range
    join <- tibble(
        `#CHROM` = "chr2",
        start = 40,
        end = 49,
        join_etc = LETTERS[1]
    )
    join_dbtbl <- db_temporary_table(db, join, "join_tbl1")
    result <- db_range_join(db, "key", "join_tbl1", "res1") |> collect()
    expect_identical(NROW(result), 1L)
    expect_identical(result$key_etc, "g")
    expect_identical(result$join_etc, "A")

    ## duplicate join range -- each range matches, so two rows
    join <- tibble(
        `#CHROM` = rep("chr2", 2),
        start = rep(40, 2),
        end = rep(49, 2),
        join_etc = LETTERS[1:2]
    )
    join_dbtbl <- db_temporary_table(db, join, "join_tbl2")
    result <- db_range_join(db, "key", "join_tbl2", "res2") |> collect()
    expect_identical(NROW(result), 2L)
    expect_true(setequal(result$key_etc, "g"))
    expect_true(setequal(result$join_etc, join$join_etc))

    db_disconnect(db)
})

test_that("'db_disconnect()' works", {
    fl <- tempfile(); dir.create(fl)
    bfc <- BiocFileCache::BiocFileCache(fl)

    db <- db_connect(ALPHAMISSENSE_RECORD, bfc)
    expect_true(DBI::dbIsValid(db))
    expect_identical(db_disconnect(db), TRUE)  # return value is invisible
    expect_identical(db_disconnect(db), FALSE) # signal db already closed
    expect_false(DBI::dbIsValid(db))
})
