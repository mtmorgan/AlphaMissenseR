#' @importFrom memoise memoise
.onLoad <-
    function(libpath, pkgname)
{
    spdl::set_level("info")
    spdl::set_pattern("* [%T][%l] %v")

    ## initialize 'ALPHAMISSENSE_STATE[["internet_available"]]
    internet_available()

    ALPHAMISSENSE_RECORD <<- Sys.getenv("ALPHAMISSENSE_RECORD", "10813168")

    am_record_json <<- memoise(am_record_json)
    am_record_is_latest <<- memoise(am_record_is_latest)
}
