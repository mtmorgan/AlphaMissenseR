#' @importFrom memoise memoise
.onLoad <-
    function(...)
{
    spdl::set_level("info")
    spdl::set_pattern("* [%T][%l] %v")

    ALPHAMISSENSE_RECORD <<- Sys.getenv("ALPHAMISSENSE_RECORD", "8360242")

    am_available <<- memoise(am_available)
}
