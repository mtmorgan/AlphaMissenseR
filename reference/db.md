# Manipulate the Database of Missense Mutations

`db_connect()` manages connections to AlphaMissense record-specific
databases. By default, connections are created once and reused.

`db_tables()` queries for the names of temporary and regular tables
defined in the database.

`db_temporary_table()` creates a temporary (for the duration of the
duckdb connection) table from a tibble.

`db_range_join()` performs a range join, finding all positions in `key`
within ranges defined by `join`. The result is stored in table `to`.

`db_disconnect()` disconnects the duckdb database and shuts down the
DuckDB server associated with the connection. Temporary tables are lost.

`db_disconnect_all()` disconnects all managed duckdb database
connection.

## Usage

``` r
db_connect(
    record = ALPHAMISSENSE_RECORD,
    bfc = BiocFileCache(),
    read_only = TRUE,
    managed = read_only
)

db_tables(db = db_connect())

db_temporary_table(db, value, to)

db_range_join(db, key, join, to)

db_disconnect(db = db_connect())

db_disconnect_all()
```

## Arguments

- record:

  character(1) Zenodo record for the AlphaMissense data resources.

- bfc:

  an object returned by `BiocFileCache()` representing the location
  where downloaded files and the parsed database will be stored. The
  default is the 'global' BiocFileCache.

- read_only:

  logical(1) open the connection 'read only'. `TRUE` protects against
  overwriting existing data and is the default.

- managed:

  logical(1) when `TRUE`, re-use an existing managed connection to the
  same database.

- db:

  `duckdb_connection` object, returned by `db_connect()`.

- value:

  a `data.frame` / `tibble` containing data to be placed in a temporary
  table, e.g., from a GenomicRanges object to be used in a range join.

- to:

  the character(1) name of the table to be created

- key:

  a character(1) table name in `db` containing missense mutation
  coordinates.

- join:

  a character(1) table name in `db` containing ranges to be used for
  joining with (filtering) `key`.

## Value

`db_connect()` returns an open `duckdb_connection` to the AlphaMissense
record-specific database.

`db_tables()` returns a character vector of database table names.

`db_temporary_table()` returns the temporary table as a dbplyr `tibble`.

`db_range_join()` returns `to` (the temporary table created from the
join) as a dbplyr tibble.

`db_disconnect()` returns `FALSE` if the connection has already been
closed or is not valid (via `dbIsValid()`) or `TRUE` if disconnection is
successful. Values are returned invisibly.

`db_disconnect_all()` returns the `db_disconnect()` value for each
connection, invisibly.

## Details

For `db_connect()`, set `managed = FALSE` when, for instance, accessing
a database in a separate process. Remember to capture the database
connection `db_unmanaged <- db_connect(managed = FALSE)` and disconnect
when done \`db_disconnect(db_unmanaged). Connections are managed by
default.

`db_temporary_table()` **overwrites** an existing table with name `to`.

`db_range_join()` **overwrites** an existing table `to`. The table `key`
is usually `"hg19"` or `"hg38"` and must have `CHROM` and `POS` columns.
The table `join` must have columns `CHROM`, `start` and `end`. Following
*Bioconductor* convention and as reported in
[`am_browse()`](https://mtmorgan.github.io/AlphaMissenseR/reference/AlphaMissense.md),
coordinates are 1-based and ranges defined by `start` and `end` are
closed. All columns from both `key` and `join` are included, so column
names (other than `CHROM`) cannot be duplicated.

`db_disconnect()` should be called on each unmanaged connection, and
once (to free the default managed connection) at the end of a session.

## Examples

``` r
db_connect()          # default 'read-only' connection
#> alphamissense_connection (read_only: TRUE; managed: TRUE)
#> record: '10813168'
#> connected: TRUE
#> tables: aa_substitutions, clinvar, gene_hg38, hg38

db_rw <- db_connect(read_only = FALSE)

am_data("hg38")       # uses the default, 'read-only' connection
#> # Source:   table<hg38> [?? x 10]
#> # Database: DuckDB 1.4.2 [mtmorgan@Darwin 25.1.0:R 4.6.0//Users/mtmorgan/Library/Caches/org.R-project.R/R/BiocFileCache/785c21cda4e0_785c21cda4e0]
#>    CHROM   POS REF   ALT   genome uniprot_id transcript_id     protein_variant
#>    <chr> <dbl> <chr> <chr> <chr>  <chr>      <chr>             <chr>          
#>  1 chr1  69094 G     T     hg38   Q8NH21     ENST00000335137.4 V2L            
#>  2 chr1  69094 G     C     hg38   Q8NH21     ENST00000335137.4 V2L            
#>  3 chr1  69094 G     A     hg38   Q8NH21     ENST00000335137.4 V2M            
#>  4 chr1  69095 T     C     hg38   Q8NH21     ENST00000335137.4 V2A            
#>  5 chr1  69095 T     A     hg38   Q8NH21     ENST00000335137.4 V2E            
#>  6 chr1  69095 T     G     hg38   Q8NH21     ENST00000335137.4 V2G            
#>  7 chr1  69097 A     G     hg38   Q8NH21     ENST00000335137.4 T3A            
#>  8 chr1  69097 A     C     hg38   Q8NH21     ENST00000335137.4 T3P            
#>  9 chr1  69097 A     T     hg38   Q8NH21     ENST00000335137.4 T3S            
#> 10 chr1  69098 C     A     hg38   Q8NH21     ENST00000335137.4 T3N            
#> # ℹ more rows
#> # ℹ 2 more variables: am_pathogenicity <dbl>, am_class <chr>
db_tables()           # connections initially share the same tables
#> [1] "aa_substitutions" "clinvar"          "gene_hg38"        "hg38"            
db_tables(db_rw)
#> [1] "aa_substitutions" "clinvar"          "gene_hg38"        "hg38"            

## ranges of interest -- the first 200000 bases on chromsomes 1-4.
ranges <- tibble(
    CHROM = paste0("chr", 1:4),
    start = rep(1, 4),
    end = rep(200000, 4)
)
db_temporary_table(db_rw, ranges, "ranges")
#> # Source:   table<ranges> [?? x 3]
#> # Database: DuckDB 1.4.2 [mtmorgan@Darwin 25.1.0:R 4.6.0//Users/mtmorgan/Library/Caches/org.R-project.R/R/BiocFileCache/785c21cda4e0_785c21cda4e0]
#>   CHROM start    end
#>   <chr> <dbl>  <dbl>
#> 1 chr1      1 200000
#> 2 chr2      1 200000
#> 3 chr3      1 200000
#> 4 chr4      1 200000

db_tables(db_rw)      # temporary table available to the db_rw connection...
#> [1] "aa_substitutions" "clinvar"          "gene_hg38"        "hg38"            
#> [5] "ranges"          
db_tables()           # ...but not to the read-only connection
#> [1] "aa_substitutions" "clinvar"          "gene_hg38"        "hg38"            

rng <- db_range_join(db_rw, "hg38", "ranges", "ranges_overlaps")
rng
#> # Source:   table<ranges_overlaps> [?? x 12]
#> # Database: DuckDB 1.4.2 [mtmorgan@Darwin 25.1.0:R 4.6.0//Users/mtmorgan/Library/Caches/org.R-project.R/R/BiocFileCache/785c21cda4e0_785c21cda4e0]
#>    CHROM   POS REF   ALT   genome uniprot_id transcript_id     protein_variant
#>    <chr> <dbl> <chr> <chr> <chr>  <chr>      <chr>             <chr>          
#>  1 chr4  59430 G     C     hg38   Q8IYB9     ENST00000610261.6 E2Q            
#>  2 chr4  59430 G     A     hg38   Q8IYB9     ENST00000610261.6 E2K            
#>  3 chr4  59431 A     C     hg38   Q8IYB9     ENST00000610261.6 E2A            
#>  4 chr4  59431 A     G     hg38   Q8IYB9     ENST00000610261.6 E2G            
#>  5 chr4  59431 A     T     hg38   Q8IYB9     ENST00000610261.6 E2V            
#>  6 chr4  59432 A     C     hg38   Q8IYB9     ENST00000610261.6 E2D            
#>  7 chr4  59432 A     T     hg38   Q8IYB9     ENST00000610261.6 E2D            
#>  8 chr4  59433 C     A     hg38   Q8IYB9     ENST00000610261.6 L3I            
#>  9 chr4  59433 C     T     hg38   Q8IYB9     ENST00000610261.6 L3F            
#> 10 chr4  59433 C     G     hg38   Q8IYB9     ENST00000610261.6 L3V            
#> # ℹ more rows
#> # ℹ 4 more variables: am_pathogenicity <dbl>, am_class <chr>, start <dbl>,
#> #   end <dbl>
rng |>
    count(CHROM) |>
    arrange(CHROM)
#> # Source:     SQL [?? x 2]
#> # Database:   DuckDB 1.4.2 [mtmorgan@Darwin 25.1.0:R 4.6.0//Users/mtmorgan/Library/Caches/org.R-project.R/R/BiocFileCache/785c21cda4e0_785c21cda4e0]
#> # Ordered by: CHROM
#>   CHROM     n
#>   <chr> <dbl>
#> 1 chr1   2018
#> 2 chr2   2028
#> 3 chr4   7488

db_disconnect(db_rw)  # explicit read-write connection
db_disconnect()       # implicit read-only connection

db_disconnect_all()
```
