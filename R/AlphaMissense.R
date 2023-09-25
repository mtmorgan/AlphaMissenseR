ALPHAMISSENSE_RECORD <- "8360242"
ALPHAMISSENSE_ZENODO <- "https://zenodo.org/"

#' @rdname AlphaMissense
#'
#' @title Retrieve AlphaMissense Resources as DuckDB Databases
#'
#' @description `am_browse()` opens a web browser at the Zenodo record
#'     for the AlphaMissense data.
#'
#' @param record character(1) Zenodo record for the AlphaMissense data
#'     resources.
#'
#' @examples
#' if (interactive())
#'     am_browse()
#'
#' @importFrom utils browseURL
#'
#' @export
am_browse <-
    function(record = ALPHAMISSENSE_RECORD)
{
    stopifnot(is_scalar_character(ALPHAMISSENSE_RECORD))

    url <- paste0(ALPHAMISSENSE_ZENODO, "/record/", record)
    browseURL(url)
}

#' @rdname AlphaMissense
#'
#' @description `am_available()` reports available datasets in the
#'     record.
#'
#' @return `am_available()` returns a tibble with columns `key`,
#'     `size`, and `link`. The meaning of key must be determined with
#'     reference to the information at `am_browse()`.
#'
#' @examples
#' am_available()
#'
#' @importFrom rjsoncons jmespath
#'
#' @importFrom dplyr tibble
#'
#' @export
am_available <-
    function(record = ALPHAMISSENSE_RECORD)
{
    url <- paste0(ALPHAMISSENSE_ZENODO, "/api/records/", record)
    json <- readLines(url, warn = FALSE)

    is_latest <- rjmespath(json, "metadata.relations.version[0].is_last")
    if (!is_latest)
        spdl::info("{} is not the most recent version", record)

    ## exclude 'README.md'
    files <- jmespath(json, "files[?type != 'md']")
    key <- rjmespath(files, "[*].key[]")
    size <- rjmespath(files, "[*].size[]")
    link <- rjmespath(files, "[*].links[].self")
    tibble(key, size, link, record)
}

#' @rdname AlphaMissense
#'
#' @description `am_data()` retrieves a single `key` from the the
#'     AlpahMissense Zenodo site and parses the file into a DuckDB
#'     database.
#'
#' @param key a character(1) 'key' from the result of
#'     `am_available()`, or a single row of the tibble returned by
#'     `am_available()`.
#'
#' @param bfc an object returned by `BiocFileCache()` representing the
#'     location where downloaded files and the parsed database will be
#'     stored. The default is the 'global' BiocFileCache.
#'
#' @param as chracter(1) type of return value.
#'
#' - `"duckdb"`: a duckdb database connection. Use
#'   `DBI::dbListTables(<duckdb_connection>)` to see available
#'   tables. Use `tbl(<duckdb_connection>, "hg28")` and similar to
#'   create a dbplyr tibble.
#'
#' - `"tbl"`: a dbplyr tbl representation of the database resource.
#'
#' - `"tsv"`: path to the tsv.gz file representing the resource and
#'   downloaded from Zenodo
#'
#' @details `am_data()` uses BiocFileCache to download and store the
#'     file and the corresponding DuckDB database.
#'
#' @return `am_data()` returns a dbplyr (database) tibble
#'     represented the downloaded and parsed file. Fields in the
#'     database are as described on the Zenodo resource page.
#'
#' @examples
#' am_data("AlphaMissense_hg38.tsv.gz")
#'
#' @importFrom dplyr tbl filter
#'
#' @importFrom rlang .data .env
#'
#' @importFrom BiocFileCache BiocFileCache bfcrpath bfcnew bfcquery
#'
#' @importFrom duckdb duckdb dbConnect dbDisconnect
#'
#' @importFrom DBI dbExecute dbListTables
#'
#' @export
am_data <-
    function(key, record = ALPHAMISSENSE_RECORD,
             bfc = BiocFileCache(),
             as = c("duckdb", "tbl", "tsv"))
{
    as <- match.arg(as)
    available <- am_available(record = record)
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
    stopifnot(NROW(key) == 1L)

    rname <- paste0("AlphaMissense_", record)
    db_tbl_name <- sub("AlphaMissense_(.*)\\.tsv\\.gz", "\\1", key$key)
    if (!NROW(bfcquery(bfc, rname)))
        ## create the BiocFileCache record
        bfcnew(bfc, rname)

    if (!NROW(bfcquery(bfc, key$link)))
        spdl::info("retrieving file '{}'", key$key)
    file_path <- bfcrpath(bfc, key$link)

    db_path <- bfcrpath(bfc, rname)
    db <- dbConnect(duckdb(db_path))
    if (!db_tbl_name %in% dbListTables(db)) {
        spdl::info("creating database table '{}'", db_tbl_name)
        sql <- paste0(
            "CREATE TABLE ", db_tbl_name, " AS ",
            "SELECT * FROM read_csv_auto('", file_path, "');"
        )
        dbExecute(db, sql)
    }
    dbDisconnect(db, shutdown = TRUE) # close to flush, etc?

    db <- dbConnect(duckdb(db_path))
    switch(
        as,
        tbl = tbl(db, db_tbl_name),
        duckdb  = db,
        tsv = bfcrpath(bfc, key$link)
    )
}
