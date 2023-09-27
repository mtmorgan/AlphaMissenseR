sql_template <-
    function(name)
{
    file <- paste0(name, ".sql")
    path <- system.file(package = "AlphaMissense", "sql", file)
    lines <- readLines(path)
    paste(lines, collapse = "\n")
}

DB_CONNECTION <- new.env(parent = emptyenv())

#' @rdname db
#'
#' @title Manipulate the Database of Missense Mutations
#'
#' @description `db_connect()` manages connections to AlphaMissense
#'     record-specific databases. By default, connections are created
#'     once and reused.
#'
#' @inheritParams am_data
#'
#' @param read_only logical(1) open the connection 'read only'.
#'     `TRUE` protects against overwriting existing data and is the
#'     default.
#'
#' @param managed logical(1) when `TRUE`, re-use an existing managed
#'     connection to the same database.
#'
#' @details For `db_connect()`, set `managed = FALSE` when, for
#'     instance, accessing a database in a separate process. Remember
#'     to capture the database connection `db_unmanaged <-
#'     db_connect(managed = FALSE)` and disconnect when done
#'     `db_disconnect(db_unmanaged). Connections opened with
#'     `read_only = FALSE` are *not* managed, and must also be
#'     disconnected manually.
#'
#' @return `db_connect()` returns an open `duckdb_connection` to the
#'     AlphaMissense record-specific database.
#'
#' @examples
#' db_connect()          # default 'read-only' connection
#'
#' db_rw <- db_connect(read_only = FALSE)
#'
#' @export
db_connect <-
    function(record = ALPHAMISSENSE_RECORD, bfc = BiocFileCache(),
             read_only = TRUE,
             managed = read_only)
{
    stopifnot(
        is_scalar_character(record),
        is_scalar_logical(managed),
        is_scalar_logical(read_only)
    )

    db_connection_name <- paste("AlphaMissense", record, read_only, sep = ":")

    create_entry <-
        ## explicitly requested...
        !managed ||
        ## ...or does not yet exist...
        !exists(db_connection_name, envir = DB_CONNECTION) ||
        ## ...or has been disconnected
        !dbIsValid(DB_CONNECTION[[db_connection_name]])
    if (create_entry) {
        spdl::info(
            "creating db connection for record '{}', read only '{}'",
            record, read_only
        )
        rname <- paste0("AlphaMissense_", record)
        if (!NROW(bfcquery(bfc, rname)))
            ## create the BiocFileCache record
            bfcnew(bfc, rname)
        db_path <- bfcrpath(bfc, rname)
        db <- dbConnect(duckdb(db_path, read_only))
        if (managed)
            DB_CONNECTION[[db_connection_name]] <- db
    } else {
        ## re-use existing connection
        db <- DB_CONNECTION[[db_connection_name]]
    }

    db
}

#' @rdname db
#'
#' @description `db_tables()` queries for the names of temporary and
#'     regular tables defined in the database.
#'
#' @param db a `duckdb_connection` object, as returned by `db_connect()`.
#'
#' @return `db_tables()` returns a character vector of database table
#'     names.
#'
#' @examples
#' am_data("hg38")       # uses the default, 'read-only' connection
#' db_tables()           # connections initially share the same tables
#' db_tables(db_rw)
#'
#' @importFrom DBI dbListTables dbWriteTable
#'
#' @export
db_tables <-
    function(db = db_connect())
{
    stopifnot(
        inherits(db, "duckdb_connection"),
        dbIsValid(db)
    )

    dbListTables(db)
}

#' @importFrom DBI dbListFields
db_table_fields <-
    function(db = db_connect(), table)
{
    stopifnot(
        inherits(db, "duckdb_connection"),
        is_scalar_character(table),
        table %in% db_tables(db)
    )

    dbListFields(db, table)
}

#' @rdname db
#'
#' @description `db_temporary_table()` creates a temporary (for the
#'     duration of the duckdb connection) table from a tibble.
#'
#' @details `db_temporary_table()` **overwrites** an existing table
#'     with name `to`.
#'
#' @param value a `data.frame` / `tibble` containing data to be placed
#'     in a temporary table, e.g., from a GenomicRanges object to be
#'     used in a range join.
#'
#' @param to the character(1) name of the table to be created
#'
#' @return `db_temporary_table()` returns the temporary table as a
#'     dbplyr `tibble`.
#'
#' @examples
#' ## ranges of interest -- the first 200000 bases on chromsomes 1-4.
#' ranges <- tibble(
#'     `#CHROM` = paste0("chr", 1:4),
#'     start = rep(1, 4),
#'     end = rep(200000, 4)
#' )
#' db_temporary_table(db_rw, ranges, "ranges")
#'
#' db_tables(db_rw)      # temporary table available to the db_rw connection...
#' db_tables()           # ...but not to the read-only connection
#'
#' @export
db_temporary_table <-
    function(db, value, to)
{
    stopifnot(
        inherits(db, "duckdb_connection"),
        is_scalar_character(to),
        inherits(value, "data.frame")
    )

    if (to %in% db_tables(db))
        spdl::info("overwriting existing table '{}'", to)

    dbWriteTable(db, to, value, temporary = TRUE, overwrite = TRUE)

    tbl(db, to)
}

#' @rdname db
#'
#' @description `db_range_join()` performs a range join, finding all
#'     positions in `key` within ranges defined by `join`. The result
#'     is stored in table `to`.
#'
#' @details
#'
#' `db_range_join()` **overwrites** and existing table with name `to`.
#' The table `key` is usually `"hg19"` or `"hg38"` and must have
#' `#CHROM` and `POS` columns. The table `join` must have columns
#' `#CHROM`, `start` and `end`. Following *Bioconductor* convention
#' and as reported in at `am_browse()`, coordinates are 1-based and
#' ranges defined by `start` and `end` are closed. All columns from
#' both `key` and `join` are included, so column names (other than
#' `#CHROM`) cannot be duplicated.
#'
#' @param key a character(1) table name in `db` containing missense
#'     mutation coordinates.
#'
#' @param join a character(1) table name in `db` containing ranges to
#'     be used for joining with (filtering) `key`.
#'
#' @return `db_range_join()` returns `to` (the temporary table created
#'     from the join) as a dbplyr tibble.
#'
#' @examples
#' rng <- db_range_join(db_rw, "hg38", "ranges", "ranges_overlaps")
#' rng
#' rng |>
#'     count(`#CHROM`) |>
#'     arrange(`#CHROM`)
#'
#' @importFrom whisker whisker.render
#'
#' @export
db_range_join <-
    function(db, key, join, to)
{
    stopifnot(
        inherits(db, "duckdb_connection"),

        is_scalar_character(key) && key %in% db_tables(db),
        all(c("#CHROM", "POS") %in% db_table_fields(db, key)),

        is_scalar_character(join) && join %in% db_tables(db),
        all(c("#CHROM", "start", "end") %in% db_table_fields(db, join)),

        is_scalar_character(to)
    )

    if (to %in% db_tables(db))
        spdl::info("overwriting existing table '{}'", to)

    spdl::info("doing range join of '{}' with '{}'", key, join)
    template <- sql_template("range_join")
    sql <- whisker.render(template)
    dbExecute(db, sql)

    tbl(db, to)
}

#' @rdname db
#'
#' @description `db_disconnect()` disconnects the duckdb database and
#'     shuts down the DuckDB server associated with the
#'     connection. Temporary tables are lost.
#'
#' @details `db_disconnect()` should be called on each unmanaged
#'     connection, and once (to free the default managed connection)
#'     at the end of a session.
#'
#' @return `db_disconnect()` returns `FALSE` if the connection has
#'     already been closed or is not valid (via `dbIsValid()`) or
#'     `TRUE` if disconnection is successful. Values are returned
#'     invisibly.
#'
#' @examples
#' db_disconnect(db_rw)  # required for each non-'managed' connection
#' db_disconnect()       # required once per session
#'
#' @importFrom DBI dbIsValid dbDisconnect
#'
#' @export
db_disconnect <-
    function(db = db_connect())
{
    stopifnot(
        inherits(db, "duckdb_connection")
    )

    result <- dbIsValid(db) && dbDisconnect(db, shutdown = TRUE)
    invisible(result)
}
