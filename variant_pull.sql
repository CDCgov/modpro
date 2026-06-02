-------------------------------------
-- COMMON TABLE EXPRESSIONS (CTEs) --
-------------------------------------

-- Extracts all variant hashes (VH; HA amino acid sequence variant) staged for HI and/or HINT antigenic testing by the "Referral System" (REFSYS) workflow 
-- LEFT ANTI JOIN removes all VHs already processed by OpenFold3 and Rosetta v3.17
-- NOTE: Listed VHs represent ONLY seasonal subtypes: "B vic", "H1 swl", and "H3"
WITH STAGED_VARIANT_HASH AS (
  SELECT
    DISTINCT
    SR.variant_hash
  FROM
    ref_sys.staged_recommendation AS SR
  LEFT ANTI JOIN
    protein_modeling.relaxed_rosetta_per_resi_energy AS RMSD_PRE ON SR.variant_hash = RMSD_PRE.variant_hash
),

-- Links staged VHs with computed "priority" score in the PROTEIN_MODELING.MODEL_PRIORITY CDP table
-- DENSE_RANK applied to re-order all extracted MODEL_PRIORITY metrics into ascending (ASC) numerical scale
-- NOTE: See MODEL_PRIORITY.SQL code for full DENSE_RANK parameters, metrics, and SQL logic
STAGED_PRIORITY AS (
  SELECT
    MP.subtype,
    MP.variant_hash,
    DENSE_RANK() OVER(
      PARTITION BY
        MP.subtype
      ORDER BY
        MP.model_priority ASC
    ) AS priority
  FROM
    STAGED_VARIANT_HASH AS SVH
  INNER JOIN
    protein_modeling.model_priority AS MP ON SVH.variant_hash = MP.variant_hash
)


----------------------
-- SELECT STATEMENT --
----------------------

SELECT
  SP_B.variant_hash
FROM
  STAGED_PRIORITY AS SP_B
WHERE
  SP_B.subtype = 'B vic'
  AND SP_B.priority <= 5
UNION
SELECT
  SP_H1.variant_hash
FROM
  STAGED_PRIORITY AS SP_H1
WHERE
  SP_H1.subtype = 'H1 swl'
  AND SP_H1.priority <= 10
UNION
SELECT
  SP_H3.variant_hash
FROM
  STAGED_PRIORITY AS SP_H3
WHERE
  SP_H3.subtype = 'H3'
  AND SP_H3.priority <= 15;
