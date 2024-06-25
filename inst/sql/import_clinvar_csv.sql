-- import csv from file at 'file_path' to data table 'db_tbl_name'
-- create the clinvar table, but with label as TEXT
CREATE TABLE {{db_tbl_name}} AS
SELECT 
    * EXCLUDE(label),
    CASE
        WHEN label = 0 THEN 'benign'
        WHEN label = 1 THEN 'pathogenic'
    END AS label
FROM read_csv_auto('{{file_path}}', delim = '{{delim}}', compression = 'gzip');
-- create levels for the label column
CREATE TYPE label_level AS ENUM ('benign', 'pathogenic');
-- coerce label to factor
ALTER TABLE {{db_tbl_name}} ALTER label TYPE label_level;
