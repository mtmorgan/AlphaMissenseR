ALPHAMISSENSE_RECORD <- "8360242"
ALPHAMISSENSE_ZENODO <- "https://zenodo.org/"

#' @rdname AlphaMissense
#'
#' @title Retrieve AlphaMissense Resources as DuckDB Databases
#'
#' @description `alphamissense_browse()` opens a web browser at the
#'     Zenodo record for the AlphaMissense data.
#'
#' @param record character(1) Zenodo record for the AlphaMissense data
#'     resources.
#'
#' @examples
#' if (interactive())
#'     alphamissense_browse()
#'
#' @importFrom utils browseURL
#'
#' @export
alphamissense_browse <-
    function(record = ALPHAMISSENSE_RECORD)
{
    stopifnot(is_scalar_character(ALPHAMISSENSE_RECORD))

    url <- paste0(ALPHAMISSENSE_ZENODO, "/record/", record)
    browseURL(url)
}

#' @rdname AlphaMissense
#'
#' @description `alphamissense_available()` reports available datasets
#'     in the record.
#'
#' @return `alphamissense_available()` returns a tibble with columns
#'     `key`, `size`, and `link`. The meaning of key must be
#'     determined with reference to the information at
#'     `alphamissense_browse()`.
#'
#' @examples
#' alphamissense_available()
#'
#' @importFrom rjsoncons jmespath
#'
#' @importFrom dplyr tibble
#'
#' @export
alphamissense_available <-
    function(record = ALPHAMISSENSE_RECORD)
{
    url <- paste0(ALPHAMISSENSE_ZENODO, "/api/records/", record)
    json <- readLines(url, warn = FALSE)

    is_latest <- identical(basename(jmespath(json, "links.latest")), record)
    if (!is_latest)
        spdl::info("{} is not the most recent version", record)

    ## exclude 'README.md'
    files <- jmespath(json, "files[?type != 'md']")
    key <- rjmespath(files, "[*].key[]")
    size <- rjmespath(files, "[*].size[]")
    link <- rjmespath(files, "[*].links[].self")
    tibble(key, size, link)
}

#' @rdname AlphaMissense
#'
#' @description `alphamissense()` retrieves a single `key` from the
#'     the AlpahMissense Zenodo site and parses the file into a DuckDB
#'     database.
#'
#' @param key a character(1) 'key' from the result of
#'     `alphamissense_available()`, or a single row of the tibble
#'     returned by `alphamissense_available()`.
#'
#' @param bfc an object returned by `BiocFileCache()` representing the
#'     location where downloaded files and the parsed database will be
#'     stored. The default is the 'global' BiocFileCache.
#'
#' @details `alphamissense()` uses BiocFileCache to store the
#'     downloaded file and the corresponding DuckDB database.
#'
#' @return `alphamissense()` returns a dbplyr (database) tibble
#'     represented the downloaded and parsed file. Fields in the
#'     database are as described on the Zenodo resource page.
#'
#' @examples
#' alphamissense("AlphaMissense_hg38.tsv.gz")
#'
#' @importFrom dplyr tbl filter
#'
#' @importFrom rlang .data .env
#'
#' @importFrom BiocFileCache BiocFileCache bfcrpath bfcnew bfcquery
#'
#' @importFrom duckdb duckdb dbConnect dbDisconnect
#'
#' @importFrom DBI dbExecute
#'
#' @export
alphamissense <-
    function(key = alphamissense_available(), bfc = BiocFileCache())
{
    available <- alphamissense_available()
    if (inherits(key, "data.frame")) {
        stopifnot(
            NROW(key) == 1L,
            all(c("key", "link") %in% colnames(key)),
            all(key$key %in% available$key)
        )
    } else {
        stopifnot(
            is.character(key),
            length(key) == 1L,
            all(key %in% available$key)
        )
        key <-
            available |>
            filter(.data$key %in% .env$key)
    }
    stopifnot(NROW(key) > 0L)

    rname <- sub(".tsv.gz", ".duckdb", key$key)
    if (!NROW(bfcquery(bfc, rname))) {
        spdl::info("retrieving file")
        path <- bfcrpath(bfc, key$link)

        spdl::info("creating database")
        db_path <- bfcnew(bfc, rname)
        db <- dbConnect(duckdb(db_path))
        sql <- paste0(
            "CREATE TABLE AlphaMissense AS ",
            "SELECT * FROM read_csv_auto('", path, "');"
        )
        dbExecute(db, sql)
        dbDisconnect(db, shutdown = TRUE)
    }
    db_path <- bfcrpath(bfc, rname)

    spdl::info("creating tbl")
    db <- dbConnect(duckdb(db_path))
    tbl(db, "AlphaMissense")
}
