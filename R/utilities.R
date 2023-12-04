ALPHAMISSENSE_STATE <- new.env(parent = emptyenv())

#' @importFrom BiocBaseUtils isCharacter isScalarCharacter
#'     isScalarLogical isScalarNumber

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
