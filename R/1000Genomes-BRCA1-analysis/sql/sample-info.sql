# Get ethnicity, gender and family relationship information for the samples.
SELECT
  Sample,
  Gender,
  Family_ID,
  Relationship,
  Population,
  Population_Description,
  Super_Population,
  Super_Population_Description
FROM
  [google.com:biggene:1000genomes.sample_info]
WHERE
  In_Phase1_Integrated_Variant_Set = TRUE
ORDER BY
  Sample
