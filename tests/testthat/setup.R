## Google Gemini suggestion to avoid duckdb message about downloaded
## extensions and secrets storage location during tests

## 1. Use a temporary subdirectory of tempdir() for all DuckDB storage
tmp_duckdb_home <- file.path(tempdir(), "duckdb_test_home")
if (!dir.exists(tmp_duckdb_home)) dir.create(tmp_duckdb_home)

## 2. Force the global R option so duckdb uses this path automatically
options(duckdb.home = tmp_duckdb_home)

## 3. Clean up the directory automatically when the test suite completes
withr::defer({
   options(duckdb.home = NULL)
   if (dir.exists(tmp_duckdb_home)) unlink(tmp_duckdb_home, recursive = TRUE)
}, teardown_env())
