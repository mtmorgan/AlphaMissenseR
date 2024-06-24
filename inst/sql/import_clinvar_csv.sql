-- import csv from file at 'file_path' to data table 'db_tbl_name'
CREATE TABLE {{db_tbl_name}} AS
SELECT *,
    CASE
        WHEN label = 0 THEN '0'
        WHEN label = 1 THEN '1'
    END AS label_category
FROM read_csv_auto('{{file_path}}', delim = '{{delim}}', compression = 'gzip');
