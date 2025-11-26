# E. Issues & Solutions

Original version: 16 October, 2023

``` r
library(AlphaMissenseR)
```

## Updating [duckdb](https://CRAN.R-project.org/package=duckdb) to 0.9.1

The R duckdb client version 0.9.1 cannot read databases created with
previous versions of the package. The duckdb error message indicates

> > am_available() Error in h(simpleError(msg, call)) : error in
> > evaluating the argument ‘drv’ in selecting a method for function
> > ‘dbConnect’: rapi_startup: Failed to open database: IO Error: Trying
> > to read a database file with version number 51, but we can only read
> > version 64.
>
> The database file was created with DuckDB version v0.8.0 or v0.8.1.
>
> The storage of DuckDB is not yet stable; newer versions of DuckDB
> cannot read old database files and vice versa. The storage will be
> stabilized when version 1.0 releases.
>
> For now, we recommend that you load the database file in a supported
> version of DuckDB, and use the EXPORT DATABASE command followed by
> IMPORT DATABASE on the current version of DuckDB.

> See the storage page for more information:
> <https://duckdb.org/internals/storage>

but in practice the most straight-forward solution is to remove existing
AlphaMissenseR data resources and ‘start again’.

The following attempts to identify AlphaMissenseR data resources cached
locally

``` r
am_rids <-
    bfcinfo() |>
    dplyr::filter(
        grepl("zenodo", rname) |
        startsWith(rname, "AlphaMissense_")
    ) |>
    pull(rid)
```

After verifying that these resources have not been created outside
AlphaMissenseR, remove them.

``` r
BiocFileCache::bfcremove(rids = am_rids)
```

Commands such as
[`am_available()`](https://mtmorgan.github.io/AlphaMissenseR/reference/AlphaMissense.md)
should report no files cached. The command

``` r
am_data("gene_hg38")
#> # Source:   table<gene_hg38> [?? x 2]
#> # Database: DuckDB 1.4.2 [mtmorgan@Darwin 25.1.0:R 4.6.0//Users/mtmorgan/Library/Caches/org.R-project.R/R/BiocFileCache/785c21cda4e0_785c21cda4e0]
#>    transcript_id      mean_am_pathogenicity
#>    <chr>                              <dbl>
#>  1 ENST00000000233.10                 0.742
#>  2 ENST00000000412.8                  0.378
#>  3 ENST00000001008.6                  0.422
#>  4 ENST00000001146.6                  0.467
#>  5 ENST00000002125.9                  0.351
#>  6 ENST00000002165.11                 0.406
#>  7 ENST00000002501.10                 0.320
#>  8 ENST00000002596.6                  0.471
#>  9 ENST00000002829.8                  0.524
#> 10 ENST00000003084.10                 0.405
#> # ℹ more rows
```

will re-download the file and insert it into a database that functions
with duckdb 0.9.1.

## Resource temporarily unavailable

Trying to access a data resource with
[`am_data()`](https://mtmorgan.github.io/AlphaMissenseR/reference/AlphaMissense.md)
may sometimes result in a DuckDB errors about “Resource unavailable”.

    > am_data("hg38")

    * [10:05:09][warning] error in evaluating the argument 'drv' in selecting a
        method for function 'dbConnect': rapi_startup: Failed to open database:
        IO Error: Could not set lock on file 
        ".../Caches/org.R-project.R/R/BiocFileCache/1ec5157ddaa2_1ec5157ddaa2":
        Resource unavailable

    Error in value[[3L]](cond) :
        failed to connect to DuckDB database, see 'Issues & Solutions' vignette

This occures when the database is being used by an independent *R*
process. The solution is to identify the process and disconnect from the
database using, e.g.,
[`db_disconnect_all()`](https://mtmorgan.github.io/AlphaMissenseR/reference/db.md).

## Finally

Remember to disconnect and shutdown all managed DuckDB connections.

``` r
db_disconnect_all()
#> * [14:48:08][info] disconnecting all registered connections
```

Database connections that are not closed correctly trigger warning
messages.

## Session information

``` r
sessionInfo()
#> R Under development (unstable) (2025-08-09 r88545)
#> Platform: aarch64-apple-darwin24.6.0
#> Running under: macOS Tahoe 26.1
#> 
#> Matrix products: default
#> BLAS:   /Users/mtmorgan/bin/R-devel/lib/libRblas.dylib 
#> LAPACK: /Users/mtmorgan/bin/R-devel/lib/libRlapack.dylib;  LAPACK version 3.12.1
#> 
#> locale:
#> [1] en_CA.UTF-8/en_CA.UTF-8/en_CA.UTF-8/C/en_CA.UTF-8/en_CA.UTF-8
#> 
#> time zone: America/Toronto
#> tzcode source: internal
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#> [1] AlphaMissenseR_1.7.1 dplyr_1.1.4         
#> 
#> loaded via a namespace (and not attached):
#>  [1] utf8_1.2.6           rappdirs_0.3.3       sass_0.4.10         
#>  [4] generics_0.1.4       spdl_0.0.5           RSQLite_2.4.4       
#>  [7] digest_0.6.39        magrittr_2.0.4       RColorBrewer_1.1-3  
#> [10] evaluate_1.0.5       grid_4.6.0           fastmap_1.2.0       
#> [13] blob_1.2.4           jsonlite_2.0.0       whisker_0.4.1       
#> [16] DBI_1.2.3            purrr_1.2.0          scales_1.4.0        
#> [19] httr2_1.2.1          textshaping_1.0.4    jquerylib_0.1.4     
#> [22] duckdb_1.4.2         cli_3.6.5            rlang_1.1.6         
#> [25] dbplyr_2.5.1         bit64_4.6.0-1        withr_3.0.2         
#> [28] cachem_1.1.0         yaml_2.3.10          BiocBaseUtils_1.12.0
#> [31] tools_4.6.0          memoise_2.0.1        ggplot2_4.0.1       
#> [34] filelock_1.0.3       curl_7.0.0           rjsoncons_1.3.2     
#> [37] vctrs_0.6.5          R6_2.6.1             BiocFileCache_3.0.0 
#> [40] lifecycle_1.0.4      RcppSpdlog_0.0.23    fs_1.6.6            
#> [43] htmlwidgets_1.6.4    bit_4.6.0            ragg_1.5.0          
#> [46] pkgconfig_2.0.3      desc_1.4.3           pkgdown_2.2.0       
#> [49] pillar_1.11.1        bslib_0.9.0          gtable_0.3.6        
#> [52] Rcpp_1.1.0           glue_1.8.0           systemfonts_1.3.1   
#> [55] xfun_0.54            tibble_3.3.0         tidyselect_1.2.1    
#> [58] dichromat_2.0-0.1    knitr_1.50           farver_2.1.2        
#> [61] htmltools_0.5.8.1    rmarkdown_2.30       compiler_4.6.0      
#> [64] S7_0.2.1
```
