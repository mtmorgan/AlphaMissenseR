test_that("am_aa_class_mode() works", {
    x <- c("ambiguous", "ambiguous", "likely_pathogenic")
    levels <- c("likely_benign", "ambiguous", "likely_pathogenic")
    expect_identical(am_aa_class_mode(x), factor("ambiguous", levels))

    x <- character()
    expect_identical(am_aa_class_mode(x), factor(x, levels))

    ## ties go to less pathogenic
    x <- c("likely_pathogenic", "ambiguous", "ambiguous", "likely_pathogenic")
    expect_identical(am_aa_class_mode(x), factor("ambiguous", levels))
})

test_that("am_aa_pathogenicity() works", {

    am_pathogenicity <- c(.1, .3, .7, .5)
    am_class <- c(
        "ambiguous", "ambiguous", "likely_pathogenic", "likely_pathogenic"
    )
    tbl <- tibble(
        uniprot_id = c(rep("P1", 3), "P2"),
        protein_variant = paste0("A1", LETTERS[2:5]),
        am_pathogenicity = am_pathogenicity,
        am_class = am_class
    )
    object <- am_aa_pathogenicity(tbl)

    pathogenicity_1 <- head(am_pathogenicity, 3L)
    pathogenicity_2 <- tail(am_pathogenicity, 1L)
    expected <- tibble(
        uniprot_id = c("P1", "P2"),
        aa_pos = rep(1L, 2),
        aa_ref = rep("A", 2),
        aa_pathogenicity_n = c(3L, 1L),
        aa_pathogenicity_mean = c(
            mean(pathogenicity_1),
            mean(pathogenicity_2)
        ),
        aa_pathogenicity_median = c(
            median(pathogenicity_1),
            median(pathogenicity_2)
        ),
        aa_pathogenicity_min = c(
            min(pathogenicity_1),
            min(pathogenicity_2)
        ),
        aa_pathogenicity_max = c(
            max(pathogenicity_1),
            max(pathogenicity_2)
        ),
        aa_pathogenicity_mode = c(
            am_aa_class_mode(head(am_class, 3L)),
            am_aa_class_mode(tail(am_class, 1L))
        )
    )

    expect_identical(object, expected)
})
