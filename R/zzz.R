#' @importFrom memoise memoise
.onLoad <-
    function(...)
{
    spdl::set_level("info")
    spdl::set_pattern("* [%T][%l] %v")

    alphamissense_available <<- memoise(alphamissense_available)
}
