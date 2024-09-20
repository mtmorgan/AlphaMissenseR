# AlphaMissenseR 1.2.0

* (v. 1.1.6) Add `gosling_plot()` for visualizing variants as bar
  or lollipop plots
* (v. 1.1.4) Add `clinvar_data()` and `clinvar_plot()` for visualizing
  ClinVar data. Merges
  <https://github.com/mtmorgan/AlphaMissenseR/pull/4> Thanks
  @tram-nguyen-n
* (v. 1.1.2) `af_predictions()` returns a tibble with 21 columns,
  instead of 20. Merges
  <https://github.com/mtmorgan/AlphaMissenseR/pull/3>.

# AlphaMissenseR 1.0.0

* (v. 0.99.21) Update Zenodo data source to record 10813168, with more
  accessible 'CC-BY-4.0' license.
* (v. 0.99.16) Housekeeping
  - Order vignettes introduction, visualization, issues.
  - Use rjsoncons (>= 1.0.1), so that jsonlite is not a hard dependency.
  - Acknowledge additional funding sources; add ImmunoOncology biocViews term.
* (v. 0.99.15) Respond to Bioconductor reviewer comments
  - Use BiocBaseUtils for input assertions.
  - Include range join SQL directly rather than via a non-exported function.
  - Report file size when downloading.
  - Improve test coverage.
  - See GitHub [issue comment][]; thanks @LiNk-NY
* (v. 0.99.11) Use an S4 class for `alphamissense_connection`,
  extending `duckdb_connection`.
* (v. 0.99.10) Introduce visualization on AlphaFold predictions.
* (v. 0.99.9) Update maintainer email address, add 'Resource unavailable' section
  to 'Issues & Solutions' vignette.
* (v. 0.99.7) Update to revised Zenodo API.
* (v. 0.99.5) BREAKING change. Requires duckdb >= 0.9.1, which cannot
  read DuckDB databases created with earlier versions. Also introduces
  changes to local cache naming. See the 'Issues & Solutions' vignette.
* (v. 0.99.4) Change back to Bioconductor 'Software' package, to allow
  for regular testing of data access code.
* (v. 0.99.1) Change to Bioconductor 'Annotation' package.
* (v. 0.99.1) Update to change in Zenodo API.

# AlphaMissenseR 0.99.0

* (v. 0.0.18) Rename (including updating existing tables) '#CHROM' to
  'CHROM' in hg19 / hg38 tables
* (v. 0.0.17) Initial Bioconductor submission.

[issue comment]: https://github.com/Bioconductor/Contributions/issues/3221#issuecomment-1804040387
