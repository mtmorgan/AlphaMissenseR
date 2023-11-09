DB_CONNECTION <- new.env(parent = emptyenv())

db_connection_id <-
    function(record, bfc, read_only)
{
    paste(
        "AlphaMissense", basename(bfccache(bfc)), record, read_only,
        sep = ":"
    )
}

db_connect_create <-
    function(record, bfc, read_only)
{
    spdl::debug(
        "creating db connection for record '{}', read only '{}'",
        record, read_only
    )
    rname <- paste0("AlphaMissense_", record)
    if (!NROW(bfcquery(bfc, rname))) {
        spdl::debug("creating new DuckDB database '{}'", rname)
        ## create the BiocFileCache record
        db_path <- bfcnew(bfc, rname)
        ## open read-write to initialize
        db0 <- dbConnect(duckdb(db_path))
        dbDisconnect(db0, shutdown = TRUE)
    }
    db_path <- bfcrpath(bfc, rname)
    db <- tryCatch({
        dbConnect(duckdb(db_path, read_only))
    }, error = function(e) {
        ## e.g., because database created by older version of duckdb
        spdl::warn("{}", conditionMessage(e))
        stop(
            "failed to connect to DuckDB database, ",
            "see 'Issues & Solutions' vignette"
        )
    })

    db
}

#' @rdname db
#'
#' @title Manipulate the Database of Missense Mutations
#'
#' @description `db_connect()` manages connections to AlphaMissense
#'     record-specific databases. By default, connections are created
#'     once and reused.
#'
#' @usage
#' db_connect(
#'     record = ALPHAMISSENSE_RECORD,
#'     bfc = BiocFileCache(),
#'     read_only = TRUE,
#'     managed = read_only
#' )
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
#' @details
#'
#' For `db_connect()`, set `managed = FALSE` when, for instance,
#' accessing a database in a separate process. Remember to capture the
#' database connection `db_unmanaged <- db_connect(managed = FALSE)`
#' and disconnect when done `db_disconnect(db_unmanaged). Connections
#' are managed by default.
#'
#' @return `db_connect()` returns an open `duckdb_connection` to the
#'     AlphaMissense record-specific database.
#'
#' @examples
#' db_connect()          # default 'read-only' connection
#'
#' db_rw <- db_connect(read_only = FALSE)
#'
#' @importFrom BiocFileCache bfccache
#'
#' @export
db_connect <-
    function(
        record = ALPHAMISSENSE_RECORD, bfc = BiocFileCache(),
        read_only = TRUE, managed = read_only)
{
    stopifnot(
        isScalarCharacter(record),
        isScalarLogical(managed),
        isScalarLogical(read_only)
    )

    id <- db_connection_id(record, bfc, read_only)

    create <-
        ## explicitly requested...
        !managed ||
        ## ...or does not yet exist...
        !exists(id, envir = DB_CONNECTION) ||
        ## ...or has been disconnected
        !dbIsValid(DB_CONNECTION[[id]])
    if (create) {
        db <- db_connect_create(record, bfc, read_only)
        if (managed)
            DB_CONNECTION[[id]] <- db
    } else {
        ## re-use existing connection
        db <- DB_CONNECTION[[id]]
    }

    alphamissense_connection(
        db,
        id = id,
        record = record,
        read_only = read_only,
        managed = managed
    )
}

db_connect_or_renew <-
    function(
        record = ALPHAMISSENSE_RECORD, bfc = BiocFileCache(),
        read_only = TRUE, managed = read_only)
{
    if (!managed) {
        return(db_connect(record, bfc, read_only, managed))
    }

    id <- db_connection_id(record, bfc, read_only)
    if (id %in% names(DB_CONNECTION)) {
        ## renew
        spdl::info(paste0(
            "renewing managed connection to record '{}',\n",
            "  read_only '{}'; previous connections are invalid"
        ), record, read_only)
        db_disconnect(db_connect(record, bfc, read_only, managed))
    }

    db_connect(record, bfc, read_only, managed)
}

#' @rdname db
#'
#' @description `db_tables()` queries for the names of temporary and
#'     regular tables defined in the database.
#'
#' @param db `duckdb_connection` object, returned by `db_connect()`.
#'
#' @return `db_tables()` returns a character vector of database table
#'     names.
#'
#' @examples
#' am_data("hg38")       # uses the default, 'read-only' connection
#' db_tables()           # connections initially share the same tables
#' db_tables(db_rw)
#'
#' @export
db_tables <-
    function(db = db_connect())
{
    stopifnot(
        is(db, "alphamissense_connection"),
        dbIsValid(db)
    )

    dbListTables(db)
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
#'     CHROM = paste0("chr", 1:4),
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
        is(db, "alphamissense_connection"),
        isScalarCharacter(to),
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
#' `db_range_join()` **overwrites** an existing table `to`.
#' The table `key` is usually `"hg19"` or `"hg38"` and must have
#' `CHROM` and `POS` columns. The table `join` must have columns
#' `CHROM`, `start` and `end`. Following *Bioconductor*
#' convention and as reported in `am_browse()`, coordinates are
#' 1-based and ranges defined by `start` and `end` are closed. All
#' columns from both `key` and `join` are included, so column names
#' (other than `CHROM`) cannot be duplicated.
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
#'     count(CHROM) |>
#'     arrange(CHROM)
#'
#' @export
db_range_join <-
    function(db, key, join, to)
{
    stopifnot(
        is(db, "alphamissense_connection"),

        isScalarCharacter(key) && key %in% db_tables(db),
        all(c("CHROM", "POS") %in% dbListFields(db, key)),

        isScalarCharacter(join) && join %in% db_tables(db),
        all(c("CHROM", "start", "end") %in% dbListFields(db, join)),

        isScalarCharacter(to)
    )

    if (to %in% db_tables(db))
        spdl::info("overwriting existing table '{}'", to)

    spdl::debug("doing range join of '{}' with '{}'", key, join)
    sql <- sql_template("range_join", key = key, join = join, to = to)
    dbExecute(db, sql)

    tbl(db, to)
}

#' @importFrom DBI dbIsValid dbDisconnect
db_disconnect_duckdb <-
    function(db)
{
    ## shut down this connection
    result <- dbIsValid(db) && dbDisconnect(db, shutdown = TRUE)
    ## remove from DB_CONNECTION environment
    is_valid <- unlist(eapply(DB_CONNECTION, dbIsValid))
    if (length(is_valid) && any(!is_valid))
        rm(list = names(is_valid)[!is_valid], envir = DB_CONNECTION)

    invisible(result)
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
#' db_disconnect(db_rw)  # explicit read-write connection
#' db_disconnect()       # implicit read-only connection
#'
#' @export
db_disconnect <-
    function(db = db_connect())
{
    stopifnot(
        is(db, "alphamissense_connection")
    )

    db_disconnect_duckdb(db)
}

#' @rdname db
#'
#' @description `db_disconnect_all()` disconnects all managed duckdb
#'     database connection.
#'
#' @return `db_disconnect_all()` returns the `db_disconnect()` value
#'     for each connection, invisibly.
#'
#' @examples
#' db_disconnect_all()
#'
#' @export
db_disconnect_all <-
    function()
{
    invisible(unlist(eapply(DB_CONNECTION, db_disconnect_duckdb)))
}
