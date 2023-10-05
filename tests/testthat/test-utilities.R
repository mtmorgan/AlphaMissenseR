test_that("'internet_available()' works", {
    state0 <- internet_available()
    on.exit(capture.output(internet_available(state0)))

    output <- capture.output({
        ## no message when state does not generate logging message
        result <- internet_available(state0)
    })
    expect_identical(result, state0)
    expect_identical(length(output), 0L)

    output <- capture.output({
        ## change of state generates logging message
        result <- internet_available(!state0)
    })
    expect_identical(result, !state0)
    expect_true(nzchar(output))
})
