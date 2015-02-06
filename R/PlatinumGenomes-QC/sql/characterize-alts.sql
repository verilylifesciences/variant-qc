# Variants are only SNPs and INDELs, with no special characters.
SELECT
  COUNT(1) AS number_of_variant_records,
  REGEXP_MATCH(alternate_bases,
    r'^[A,C,G,T]+$') AS alt_contains_no_special_characters,
  MAX(LENGTH(reference_bases)) AS max_ref_len,
  MAX(LENGTH(alternate_bases)) AS max_alt_len
FROM
  [_THE_TABLE_]
OMIT RECORD IF EVERY(alternate_bases IS NULL)
GROUP BY
  alt_contains_no_special_characters
