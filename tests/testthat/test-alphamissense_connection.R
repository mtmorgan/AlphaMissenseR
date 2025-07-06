test_that("show(<alphamissense_connection>) works", {
    db <- db_connect()
    expect_output(show(db), "read_only: TRUE; managed: TRUE")
    expect_output(show(db), "connected: TRUE")
    expect_output(show(db), "tables: ")

    db <- db_connect(read_only = FALSE, managed = FALSE)
    expect_output(show(db), "read_only: FALSE; managed: FALSE")
    db_disconnect(db)
    expect_output(show(db), "connected: FALSE")
    expect_output(show(db), "tables: $")

    db_disconnect_all()
})

test_that("duckdb_connection is unchanged from duckdb 1.3.1; see https://github.com/mtmorgan/AlphaMissenseR/issues/12", {
    expect_setequal(
        names(getSlots("duckdb_connection")),
        c("conn_ref", "driver", "debug", "convert_opts", "reserved_words", 
            "timezone_out", "tz_out_convert", "bigint"
        )
    )
    expect_setequal(
        names(getSlots("alphamissense_connection")),
        c("id", "record", "read_only", "managed", "conn_ref", "driver", 
            "debug", "convert_opts", "reserved_words", "timezone_out",
            "tz_out_convert", "bigint"
        )
    )
})
