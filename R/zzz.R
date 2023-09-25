#' @importFrom memoise memoise
.onLoad <-
    function(...)
{
    spdl::set_level("info")
    spdl::set_pattern("* [%T][%l] %v")

    am_available <<- memoise(am_available)
}
