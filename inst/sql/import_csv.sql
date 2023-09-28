-- import csv from file at 'file_path' to data table 'db_tbl_name'
CREATE TABLE {{db_tbl_name}} AS
SELECT * FROM read_csv_auto('{{file_path}}');
