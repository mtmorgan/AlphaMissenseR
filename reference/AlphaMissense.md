# Retrieve AlphaMissense Resources as DuckDB Databases

`ALPHAMISSENSE_RECORD` is a constant identifier corresponding to the
default version of the AlphaMissense resource to use.

`am_browse()` opens a web browser at the Zenodo record for the
AlphaMissense data.

`am_available()` reports available datasets in the record.

`am_data()` retrieves a single `key` from the the AlpahMissense Zenodo
site and parses the file into a DuckDB database.

## Usage

``` r
ALPHAMISSENSE_RECORD

am_browse(record = ALPHAMISSENSE_RECORD)

am_available(record = ALPHAMISSENSE_RECORD, bfc = BiocFileCache())

am_data(
    key,
    record = ALPHAMISSENSE_RECORD,
    bfc = BiocFileCache(),
    as = c("tbl", "tsv")
)
```

## Format

An object of class `character` of length 1.

## Arguments

- record:

  character(1) Zenodo record for the AlphaMissense data resources.

- bfc:

  an object returned by `BiocFileCache()` representing the location
  where downloaded files and the parsed database will be stored. The
  default is the 'global' BiocFileCache.

- key:

  a character(1) 'key' from the result of `am_available()`, or a single
  row of the tibble returned by `am_available()`.

- as:

  chracter(1) type of return value.

  - `"tbl"`: a dbplyr tbl representation of the database resource.

  - `"tsv"`: path to the tsv.gz file representing the resource and
    downloaded from Zenodo

## Value

`am_available()` returns a tibble with columns `key`, `size`, and
`link`. The meaning of key must be determined with reference to the
information at `am_browse()`.

`am_data()` returns a dbplyr (database) tibble represented the
downloaded and parsed file. Fields in the database are as described on
the Zenodo resource page.

## Details

`ALPHAMISSENSE_RECORD` can be set **before the package is loaded** with
the environment variable of the same name, e.g.,
`Sys.setenv(ALPHAMISSENSE_RECORD = "10813168")`. The default is the most
recent version (version 3) as checked on 11 April, 2024.

`am_data()` uses BiocFileCache to download and store the file and the
corresponding DuckDB database.

## Examples

``` r
ALPHAMISSENSE_RECORD
#> [1] "10813168"

if (interactive())
    am_browse()

am_available()
#> # A tibble: 7 × 6
#>   record   key                             size cached filename            link 
#>   <chr>    <chr>                          <dbl> <lgl>  <chr>               <chr>
#> 1 10813168 gene_hg38                     253636 TRUE   AlphaMissense_gene… http…
#> 2 10813168 isoforms_hg38             1177361934 FALSE  AlphaMissense_isof… http…
#> 3 10813168 isoforms_aa_substitutions 2461351945 FALSE  AlphaMissense_isof… http…
#> 4 10813168 hg38                       642961469 TRUE   AlphaMissense_hg38… http…
#> 5 10813168 hg19                       622293310 FALSE  AlphaMissense_hg19… http…
#> 6 10813168 gene_hg19                     243943 FALSE  AlphaMissense_gene… http…
#> 7 10813168 aa_substitutions          1207278510 TRUE   AlphaMissense_aa_s… http…

am_data("hg38")
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

## close the connection opened when adding the data
db_disconnect()
```
