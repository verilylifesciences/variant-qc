#standardSQL
--
-- Check whether variants are only SNPs and INDELs, with no special characters.
--
SELECT
  COUNT(1) AS number_of_variant_records,
  REGEXP_CONTAINS(alt,
             r'^[ACGT]+$') AS alt_contains_no_special_characters,
  MAX(LENGTH(reference_bases)) AS max_ref_len,
  MAX(LENGTH(alt)) AS max_alt_len
FROM
  `@GENOME_CALL_OR_MULTISAMPLE_VARIANT_TABLE` v, v.alternate_bases alt
WHERE
  # Skip non-variant segments.
  EXISTS (SELECT alt FROM UNNEST(v.alternate_bases) alt WHERE alt NOT IN ("<NON_REF>", "<*>"))
GROUP BY
  alt_contains_no_special_characters
