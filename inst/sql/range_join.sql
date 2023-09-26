-- range join of 'key' with 'join'; overwrite any existing table 'to'
DROP TABLE IF EXISTS {{to}};
CREATE TEMP TABLE {{to}} AS
SELECT
    {{key}}.*,
    {{join}}.* EXCLUDE ('#CHROM')
FROM {{key}}
JOIN {{join}}
    ON {{join}}."#CHROM" = {{key}}."#CHROM"
    AND {{join}}.start <= {{key}}.POS
    AND {{join}}.end >= {{key}}.POS;
