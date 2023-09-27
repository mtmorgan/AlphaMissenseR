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
