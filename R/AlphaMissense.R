#' @rdname AlphaMissense
#'
#' @title Retrieve AlphaMissense Resources as DuckDB Databases
#'
#' @description `ALPHAMISSENSE_RECORD` is a constant identifier
#'     corresponding to the default version of the AlphaMissense
#'     resource to use.
#'
#' @details `ALPHAMISSENSE_RECORD` can be set **before the package is
#'     loaded** with the environment variable of the same name, e.g.,
#'     `Sys.setenv(ALPHAMISSENSE_RECORD = "8208688")`. The default is
#'     the most recent version (version 2) on 25 September, 2023.
#'
#' @examples
#' ALPHAMISSENSE_RECORD
#'
#' @export
ALPHAMISSENSE_RECORD <- NULL # set in .onLoad()

ALPHAMISSENSE_ZENODO <- "https://zenodo.org"

#' @rdname AlphaMissense
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
    stopifnot(isScalarCharacter(record))

    url <- paste0(ALPHAMISSENSE_ZENODO, "/record/", record)
    browseURL(url)
}

## memoized in .onLoad
am_record_json <-
    function(record)
{
    url <- paste0(ALPHAMISSENSE_ZENODO, "/api/records/", record)
    readLines(url, warn = FALSE)
}

am_data_license <-
    function(record)
{
    if (internet_available()) {
        json <- am_record_json(record)
        license <- jmespath(json, "metadata.license.id")
    } else {
        license <- "unknown (internet not available)"
    }
    spdl::info("data licensed under '{}'", toupper(license))
    license
}

am_record_is_latest <-
    function(record)
{
    json <- am_record_json(record)
    latest_version_url <- jmespath(json, "links.latest")
    latest_version <- readLines(latest_version_url, warn = FALSE)
    latest_version_id <- jmespath(json, "id")
    identical(record, latest_version_id)
}

am_available_from_internet <-
    function(record, bfc)
{
    spdl::debug("am_available_from_internet()")
    json <- am_record_json(record)
    if (!am_record_is_latest(record))
        spdl::info("{} is not the most recent version", record)

    ## exclude 'README.md'
    files <- jmespath(json, "files[?!ends_with(key, '.md')]")
    filename <- rjmespath(files, "[*].key")
    key <- sub("AlphaMissense_(.*)\\.tsv\\.gz", "\\1", filename)
    size <- rjmespath(files, "[*].size")
    db <- db_connect(record, bfc, managed = FALSE)
    cached <- key %in% db_tables(db)
    db_disconnect(db)
    link <- rjmespath(files, "[*].links[].self")
    link <- sub("(.*/files/).*", "\\1", link)
    link <- paste0(link, filename, "/content")
    tibble(record, key, size, cached, filename, link)
}

#' @importFrom BiocFileCache bfcinfo
am_available_from_cache <-
    function(record, bfc)
{
    spdl::debug("am_available_from_cache()")
    query <- paste0("AlphaMissense_", record)
    tbl <- bfcquery(query = query)
    if (!NROW(tbl))
        stop("no off-line resources available")
    db <- db_connect(record, bfc, managed = FALSE)
    key <- db_tables(db)
    db_disconnect(db)
    n_rows <- length(key)
    tibble(
        record = rep(record, n_rows),
        key,
        size = rep(NA_integer_, n_rows),
        cached = rep(TRUE, n_rows),
        filename = paste0("AlphaMissense_", key, ".tsv.gz"),
        link = rep(NA_character_, n_rows)
    )
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
#' @importFrom dplyr tibble mutate select .data
#'
#' @export
am_available <-
    function(record = ALPHAMISSENSE_RECORD, bfc = BiocFileCache())
{
    stopifnot(isScalarCharacter(record))

    if (internet_available()) {
        am_available_from_internet(record, bfc)
    } else {
        am_available_from_cache(record, bfc)
    }
}

am_data_import_csv <-
    function(record, bfc, db_tbl_name, rname, fpath)
{
    ## must be 'read_only = FALSE' so new database can be created. use
    ## 'managed = FALSE' so connection is independent of other
    ## read-write connections
    db_rw <- db_connect(record, bfc, read_only = FALSE)
    renew <- FALSE
    if (!db_tbl_name %in% db_tables(db_rw)) {
        spdl::info("downloading or finding local file")
        file_path <- bfcrpath(bfc, rnames  = rname, fpath = fpath)
        spdl::info("creating database table '{}'", db_tbl_name)
        sql <- sql_template(
            "import_csv", db_tbl_name = db_tbl_name, file_path = file_path
        )
        dbExecute(db_rw, sql)

        renew <- TRUE
    }
    if ("#CHROM" %in% dbListFields(db_rw, db_tbl_name)) {
        spdl::info("renaming '#CHROM' to 'CHROM' in table '{}'", db_tbl_name)
        sql <- sql_template(
            "rename_column",
            db_tbl_name = db_tbl_name, from = "#CHROM", to = "CHROM"
        )
        dbExecute(db_rw, sql)

        renew <- TRUE
    }
    db_disconnect(db_rw)

    if (renew) {
        ## flush managed read-only connection
        ## FIXME: but this invalidates existing read-only connections
        db <- db_connect_or_renew(record, bfc)
    } else {
        db <- db_connect(record, bfc)
    }
    tbl(db, db_tbl_name)
}

#' @rdname AlphaMissense
#'
#' @description `am_data()` retrieves a single `key` from the the
#'     AlpahMissense Zenodo site and parses the file into a DuckDB
#'     database.
#'
#' @usage
#' am_data(
#'     key,
#'     record = ALPHAMISSENSE_RECORD,
#'     bfc = BiocFileCache(),
#'     as = c("tbl", "tsv")
#' )
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
#' am_data("hg38")
#'
#' ## close the connection opened when adding the data
#' db_disconnect()
#'
#' @importFrom dplyr tbl filter
#'
#' @importFrom rlang .data .env
#'
#' @importFrom BiocFileCache BiocFileCache bfcrpath bfcnew bfcquery
#'
#' @importFrom duckdb duckdb dbConnect dbDisconnect
#'
#' @export
am_data <-
    function(
        key,
        record = ALPHAMISSENSE_RECORD, bfc = BiocFileCache(),
        as = c("tbl", "tsv"))
{
    as <- match.arg(as)
    available <- am_available(record = record)
    if (inherits(key, "data.frame")) {
        stopifnot(
            NROW(key) == 1L,
            all(c("key", "filename", "link") %in% colnames(key)),
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

    if (!NROW(bfcquery(bfc, key$filename))) {
        size <- structure(key$size, class = "object_size")
        spdl::info(
            "retrieving file name '{}' ({})",
            key$filename, format(size, units = "auto")
        )
        am_data_license(record)
    }

    db_rname <- paste0("AlphaMissense_", record)
    db_tbl_name <- key$key
    if (!NROW(bfcquery(bfc, db_rname)))
        ## create the BiocFileCache record
        bfcnew(bfc, db_rname)
    tbl <- am_data_import_csv(record, bfc, db_tbl_name, key$filename, key$link)

    switch(
        as,
        tbl = tbl,
        tsv = bfcrpath(bfc, key$link)
    )
}
