#' @importFrom methods setClass new is
alphamissense_connection <- setClass(
    "alphamissense_connection",
    contains = "duckdb_connection",
    slots = c(
        id = "character", record = "character",
        read_only = "logical", managed = "logical"
    )
)

## non-public accessors
record <- function(x) x@record

read_only <- function(x) x@read_only

managed <- function(x) x@managed

#' @importFrom DBI dbIsValid dbListTables dbListFields dbExecute
#'     dbWriteTable
setMethod("show", "alphamissense_connection", function(object) {
    status <- paste0(
        "read_only: ", read_only(object), "; ",
        "managed: ", managed(object)
    )

    is_valid <- dbIsValid(object)

    tables <- character()
    if (is_valid) {
        tables <- toString(dbListTables(object))
    }

    cat(
        class(object), " (", status, ")\n",
        "record: '", record(object), "'\n",
        "connected: ", is_valid, "\n",
        "tables: ", tables, "\n",
        sep = ""
    )
})
