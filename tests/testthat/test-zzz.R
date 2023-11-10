test_that(".onLoad() initializes variables and memoises functions", {
    expect_true(nzchar(ALPHAMISSENSE_RECORD))
    expect_true(isScalarLogical(ALPHAMISSENSE_STATE[["internet_available"]]))
    expect_true(memoise::is.memoised(am_record_json))
    expect_true(memoise::is.memoised(am_record_is_latest))
})
