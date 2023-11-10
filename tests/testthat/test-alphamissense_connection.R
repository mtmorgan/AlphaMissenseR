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
