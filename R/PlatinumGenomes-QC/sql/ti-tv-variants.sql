#Show the location of each mutation
SELECT
      reference_name,
      start,
      reference_bases,
      alternate_bases
    FROM
      [genomics-public-data:platinum_genomes.variants]
    WHERE
      reference_name = 'chr17'
      AND start BETWEEN 41196311
      AND 41277499     
      AND LENGTH(alternate_bases) == 1
      AND LENGTH(reference_bases) == 1
    GROUP BY
      reference_name,
      start,
      end,
      reference_bases,
      alternate_bases
    ORDER BY
      start