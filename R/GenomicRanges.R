#' @rdname GenomicRanges
#'
#' @title Create GenomicRanges Objects from AlphaMissense Annotations
#'
#' @description `to_GPos()` coerces a tibble derived from
#'     `am_data("hg19")` or `am_data("hg38")` resources to
#'     GenomicRanges `GPos` objects.
#'
#' @param tbl a tibble derived from `am_data("hg19")` or
#'     `am_data("hg38")`. The tibble must have columns `CHROM`, `POS`,
#'     and `genome`.
#'
#' @return `to_GPos()` returns a `GPos` object, which can be used in
#'     the same was a `GRanges` object for range-based filtering and
#'     annotation
#'
#' @examples
#' am_data("hg38") |>
#'     filter(CHROM == "chr2", POS < 10000000, REF == "G") |>
#'     select(-REF) |>
#'     to_GPos()
#'
#' @importFrom dplyr distinct collect pull
#'
#' @export
to_GPos <-
    function(tbl)
{
    stopifnot(
        inherits(tbl, "tbl") || inherits(tbl, "data.frame"),
        all(c("CHROM", "POS", "genome") %in% colnames(tbl))
    )

    if (!requireNamespace("GenomicRanges", quietly = TRUE))
        stop(
            "use 'BiocManager::install(\"GenomicRanges\")' to support ",
            "'to_GPos()'"
        )

    genome <-
        tbl |>
        distinct(.data$genome) |>
        pull()
    if (length(genome) > 1L)
        stop("'to_GPos()' 'genome' field must have exactly one distinct value")
    seqinfo <- tryCatch({
        GenomeInfoDb::Seqinfo(genome = genome) |>
            GenomeInfoDb::keepStandardChromosomes()
    }, error = function(e) {
        spdl::warn(
            "'to_GPos()' failed to discover 'seqinfo': {}",
            conditionMessage(e)
        )
        NULL
    })

    tbl <- collect(tbl)
    column_names <- colnames(tbl)

    gpos_args <- c(
        list(
            seqnames = pull(tbl, "CHROM"),
            pos = pull(tbl, "POS"),
            seqinfo = seqinfo
        ), 
        tbl |> select(-c("CHROM", "POS", "genome")) |> as.list()
    )
    do.call(GenomicRanges::GPos, gpos_args)
}
