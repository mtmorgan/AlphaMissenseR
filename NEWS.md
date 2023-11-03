# AlphaMissenseR 1.0.0

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

# AlphaMissense 0.99.0

* (v. 0.0.18) Rename (including updating existing tables) '#CHROM' to
  'CHROM' in hg19 / hg38 tables
* (v. 0.0.17) Initial Bioconductor submission.
