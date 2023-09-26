sql_template <-
    function(name)
{
    file <- paste0(name, ".sql")
    path <- system.file(package = "AlphaMissense", "sql", file)
    lines <- readLines(path)
    paste(lines, collapse = "\n")
}

#' @rdname db
#'
#' @title Manipulate the Database of Missense Mutations
#'
#' @description `db_tables()` queries for the names of temporary and
#'     regular tables defined in the database.
#'
#' @param db a `duckdb_connection` object, as returned by `am_data()`.
#'
#' @return `db_tables()` returns a character vector of database table
#'     names.
#'
#' @examples
#' db <- am_data("hg38")
#' db_tables(db)
#'
#' @importFrom DBI dbListTables dbWriteTable
#'
#' @export
db_tables <-
    function(db)
{
    stopifnot(
        inherits(db, "duckdb_connection"),
        dbIsValid(db)
    )

    dbListTables(db)
}

#' @importFrom DBI dbListFields
db_table_fields <-
    function(db, table)
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
#' db_temporary_table(db, ranges, "ranges")
#'
#' db_tables(db)
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
#' `db_range_join()` **overwites** and existing table with name `to`.
#'
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
#' rng <- db_range_join(db, "hg38", "ranges", "ranges_overlaps")
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
#' @return `db_disconnect()` returns `FALSE` if the connection has
#'     already been closed or is not valid (via `dbIsValid()`) or
#'     `TRUE` if disconnection is successful. Values are returned
#'     invisibly.
#'
#' @examples
#' db_disconnect(db)
#'
#' @importFrom DBI dbIsValid dbDisconnect
#'
#' @export
db_disconnect <-
    function(db)
{
    stopifnot(
        inherits(db, "duckdb_connection")
    )

    result <- dbIsValid(db) && dbDisconnect(db, shutdown = TRUE)
    invisible(result)
}
