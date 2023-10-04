test_that("to_GPos() works", {
    genome <- GenomeInfoDb::genome
    seqinfo <- GenomicRanges::seqinfo
    seqnames <- GenomicRanges::seqnames
    pos <- GenomicRanges::pos
    mcols <- GenomicRanges::mcols
    NROW <- BiocGenerics::NROW

    chr <- paste0("chr", rep(1:2, c(3, 5)))
    tbl <- tibble(
        `#CHROM` = chr,
        POS = 1:8,
        extra = letters[1:8],
        genome = "hg38"
    )
    expect_no_condition(gpos <- to_GPos(tbl))
    expect_identical(NROW(gpos), 8L)
    expect_identical(as.character(seqnames(gpos)), chr)
    expect_identical(pos(gpos), 1:8)
    expect_identical(unique(genome(gpos)), "hg38")
    expect_identical(names(mcols(gpos)), "extra")
    expect_identical(gpos$extra, letters[1:8])
    expect_identical(names(genome(gpos)), paste0("chr", c(1:22, "X", "Y", "M")))

    tbl <- tibble(`#CHROM` = character(), POS = integer(), genome = character())
    expect_no_condition(gpos <- to_GPos(tbl))
    expect_identical(NROW(gpos), 0L)
    expect_identical(genome(gpos), stats::setNames(nm = character()))
    expect_identical(colnames(mcols(gpos)), character())

    tbl <- tibble(
        `#CHROM` = character(), POS = integer(), genome = character(),
        extra = character()
    )
    expect_no_condition(gpos <- to_GPos(tbl))
    expect_identical(colnames(mcols(gpos)), "extra")
    expect_identical(mcols(gpos)$extra, character(0))

    tbl <- tibble(
        `CHROM` = character(), POS = integer(), genome = character()
    )
    expect_error(to_GPos(tbl))

    tbl <- tibble(
        `#CHROM` = character(), pos = integer(), genome = character()
    )
    expect_error(to_GPos(tbl))
    tbl <- tibble()
    expect_error(to_GPos(tbl))
})
