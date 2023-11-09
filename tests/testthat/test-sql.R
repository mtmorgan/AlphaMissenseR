test_that("sql templates exist with expected variables", {
    ## 'range_join' args 'key', 'join', 'to'
    sql <- sql_template("range_join")
    expect_true(isScalarCharacter(sql))

    sql <- sql_template(
        "range_join",
        key = "key", join = "join", to = "to",
        red_herring = "red_herring"
    )
    expect_true(
        grepl("key", sql) &&
        grepl("join", sql) &&
        grepl("to", sql) &&
        !grepl("red_herring", sql)
    )

    ## 'import_csv' args 'db_tbl_name', 'file_path'
    sql <- sql_template("import_csv")
    expect_true(isScalarCharacter(sql))

    sql <- sql_template(
        "import_csv",
        db_tbl_name = "db_tbl_name", file_path = "file_path", 
        red_herring = "red_herring"
    )
    expect_true(
        grepl("db_tbl_name", sql) &&
        grepl("file_path", sql) &&
        !grepl("red_herring", sql)
    )
})
