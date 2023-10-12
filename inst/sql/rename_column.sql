-- rename db_tbl_name column 'from' to 'to'
-- primary use is to rename "#CHROM" to "CHROM"; 'from' needs to be quoted
ALTER TABLE {{db_tbl_name}}
RENAME COLUMN "{{from}}" TO {{to}};
