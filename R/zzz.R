#' @importFrom memoise memoise
.onLoad <-
    function(libpath, pkgname)
{
    spdl::set_level("info")
    spdl::set_pattern("* [%T][%l] %v")

    ALPHAMISSENSE_RECORD <<- Sys.getenv("ALPHAMISSENSE_RECORD", "8360242")

    am_record_json <<- memoise(am_record_json)
}
