#standardSQL
--
-- Check whether variants are only SNPs and INDELs, with no special characters.
--
SELECT
  COUNT(1) AS number_of_variant_records,
  REGEXP_CONTAINS(a.alt,
             r'^[ACGT]+$') AS alt_contains_no_special_characters,
  MAX(LENGTH(reference_bases)) AS max_ref_len,
  MAX(LENGTH(a.alt)) AS max_alt_len
FROM
  `{{ GENOME_CALL_OR_MULTISAMPLE_VARIANT_TABLE }}` AS v, UNNEST(v.alternate_bases) AS a
WHERE
  # Skip alternate bases for non-variant segments or "any" alternate.
  a.alt NOT IN ("<NON_REF>", "<*>")
GROUP BY
  alt_contains_no_special_characters
