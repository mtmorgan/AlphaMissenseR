ALPHAMISSENSE_STATE <- new.env(parent = emptyenv())

is_scalar_character <-
    function(x)
{
    is.character(x) && length(x) == 1L && !is.na(x)
}

is_scalar_logical <-
    function(x)
{
    is.logical(x) && length(x) == 1L && !is.na(x)
}

#' @importFrom rjsoncons jmespath
#'
#' @importFrom jsonlite parse_json
rjmespath <-
    function(data, path, ..., simplifyVector = TRUE)
{
    json <- jmespath(data, path, ...)
    parse_json(json, simplifyVector = simplifyVector)
}

#' @importFrom curl has_internet
internet_available <-
    function(available = has_internet())
{
    ## use ALPHAMISSENSE_STATE to only print info on change-of-state
    state <- ALPHAMISSENSE_STATE[["internet_available"]]
    if (!identical(available, state)) {
        if (!is.null(state))
            spdl::info("internet available: {}", available)
        ALPHAMISSENSE_STATE[["internet_available"]] <- available
    }
    ALPHAMISSENSE_STATE[["internet_available"]]
}