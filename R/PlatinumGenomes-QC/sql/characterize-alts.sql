# Check whether variants are only SNPs and INDELs, with no special characters.
SELECT
  COUNT(1) AS number_of_variant_records,
  REGEXP_MATCH(alternate_bases,
    r'^[A,C,G,T]+$') AS alt_contains_no_special_characters,
  MAX(LENGTH(reference_bases)) AS max_ref_len,
  MAX(LENGTH(alternate_bases)) AS max_alt_len
FROM
  [_THE_TABLE_]
# In some datasets, alternate_bases will be empty (therefore NULL) for non-variant segments.
# In other datasets, alternate_bases will have the value "<NON_REF>" for non-variant segments.
OMIT RECORD IF EVERY(alternate_bases IS NULL) OR EVERY(alternate_bases = "<NON_REF>")
GROUP BY
  alt_contains_no_special_characters
