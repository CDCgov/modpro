------------------------------------
-- COMMON TABLE EXPRESSIONS (CTEs) --
-------------------------------------

-- Extract all unique HA variant sequences archived in DAIS and PUBSEQ
-- See TAXONOMY.SQL query in "Referral System (REFSYS)" pipeline for more details on applied contextual filters
WITH LINK_AA_SEQUENCE AS (
  SELECT
    TAXA.ctype,
    TAXA.subtype,
    TAXA.genetic_group,
    TAXA.clade,
    TAXA.protein,
    TAXA.variant_hash,
    TAXA.rep_nt_id,
    TAXA.cds_seq_length,
    TAXA.aa_seq_length,
    TAXA.shift_indel,
    TAXA.has_insertion,
    ALN.aa_seq,
    TAXA.variant_substitutions,
    TAXA.variant_base_substitutions
  FROM
    ref_sys.taxonomy AS TAXA
  INNER JOIN
    dais.alignments AS ALN ON TAXA.variant_hash = ALN.variant_hash
    AND TAXA.rep_nt_id = ALN.nt_id
  WHERE
    (
      TAXA.consensus_variants IS NULL
      OR TAXA.consensus_variants = 0
    )
),

-- Link specimen COUNTs for both internal (i.e., DISC) and external (i.e., PUBSEQ/GISAID) surveillance
-- Contexts include the following:
--   (1) CUMULATIVE - Includes both internal/external specimen COUNTs with no temporal restrictions
--   (2) WINDOWED - Internal/external specimen COUNTs observed within the last 274 days (nine months)
LINK_SURVEILLANCE_METRICS AS (
  SELECT
    LAAS.ctype,
    LAAS.subtype,
    LAAS.genetic_group,
    LAAS.clade,
    LAAS.protein,
    LAAS.variant_hash,
    LAAS.rep_nt_id,
    SURV.cumulative_surveillance,
    SURV.windowed_surveillance,
    SURV.internal_surveillance,
    SURV.external_surveillance,
    SURV.last_seen,
    LAAS.cds_seq_length,
    LAAS.aa_seq_length,
    LAAS.shift_indel,
    LAAS.has_insertion,
    LAAS.aa_seq,
    LAAS.variant_substitutions,
    LAAS.variant_base_substitutions
  FROM
    LINK_AA_SEQUENCE AS LAAS
  INNER JOIN
    ref_sys.surveillance AS SURV ON LAAS.variant_hash = SURV.variant_hash
),

-- Enumerate number of amino acid substitutions in relevant antigenic sites, binding domains, etc.
DETECT_SIGNATURE_SEQUENCE_VARIANTS AS (
  SELECT
    DFEAT.variant_hash,
    SUM(DFEAT.element_detected) AS signature_changes
  FROM (
    SELECT
      DISTINCT
      LSM.variant_hash,
      FBS.site,
      udx.contains_element(
        regexp_replace(
          LSM.variant_base_substitutions,
          '[A-Z\-]',
          ''
        ),
        CAST(FBS.site AS STRING),
        ','
      ) AS element_detected
    FROM
      LINK_SURVEILLANCE_METRICS AS LSM
    INNER JOIN
      protein_modeling.feature_by_site AS FBS ON LSM.protein = FBS.protein
      AND IF(
        (
          LSM.subtype = 'B vic'
          AND LSM.protein = 'HA'
        ),
        FBS.reference_id = 'BRISBANE60',
        FBS.reference_id != 'BRISBANE60'
      )
      AND FBS.feature_type IN (
        'binding',
        'modification'
      )
      AND FBS.feature IN (
        '120',
        '150',
        '160',
        '190',
        'Ca',
        'Cb',
        'Sa',
        'Sb',
        'antigenic_site_A',
        'antigenic_site_B',
        'antigenic_site_C',
        'antigenic_site_D',
        'antigenic_site_E',
        'receptor_binding_site',
        'n_linked_glycosylation_site',
        'koeln_major'
      )
    ) AS DFEAT
    GROUP BY
      DFEAT.variant_hash
),

-- Initial aggregation of "valid", seasonal HI protocol logfold antigenic records at the antigen specimen-level (i.e., AG_CDC_ID)
-- See FILTERED_ANTIGENIC_RESPONSES.SQL query in "Referral System (REFSYS)" pipeline for more details on applied contextual filters
LINK_HI_ANTIGENIC_RECORD AS (
  SELECT
    DISTINCT
    LAAS.subtype,
    ANTS.test_protocol,
    ANTS.sr_variant_hash,
    ANTS.lot,
    LAAS.variant_hash,
    ANTS.n AS testedness,
    ANTS.n_test_id,
    ANTS.avg_logfold
  FROM
    LINK_AA_SEQUENCE AS LAAS
  INNER JOIN
    ref_sys.antigenic_summaries AS ANTS ON LAAS.subtype = ANTS.subtype
    AND LAAS.variant_hash = ANTS.ag_variant_hash
  WHERE
    ANTS.aggregation = 'reference'
    AND ANTS.test_protocol IN (
      'hi_protocol',
      'hi_oseltamivir_protocol'
    )
    AND NOT ANTS.egg_passage
),

-- Initial aggregation of "valid", seasonal HINT protocol logfold antigenic records at the antigen specimen-level (i.e., AG_CDC_ID)
-- See FILTERED_ANTIGENIC_RESPONSES.SQL query in "Referral System (REFSYS)" pipeline for more details on applied contextual filters
LINK_HINT_ANTIGENIC_RECORD AS (
  SELECT
    DISTINCT
    LAAS.subtype,
    ANTS.test_protocol,
    ANTS.sr_variant_hash,
    ANTS.lot,
    LAAS.variant_hash,
    ANTS.n AS testedness,
    ANTS.n_test_id,
    ANTS.avg_logfold
  FROM
    LINK_AA_SEQUENCE AS LAAS
  INNER JOIN
    ref_sys.antigenic_summaries AS ANTS ON LAAS.subtype = ANTS.subtype
    AND LAAS.variant_hash = ANTS.ag_variant_hash
  WHERE
    ANTS.aggregation = 'reference'
    AND ANTS.test_protocol = 'hint_protocol'
    AND NOT ANTS.egg_passage
),

-- Apply categorical "bins" to average, logfold HI/HINT responses: LIKE, LIKE(ISH), and LOW (see below for criteria)
BIN_ANTIGENIC_TRENDS AS (
  SELECT
    DISTINCT
    LAAS.subtype,
    LAAS.protein,
    LAAS.variant_hash,
    LAAS.rep_nt_id,
    HIAR.testedness AS hi_testedness,
    HIAR.n_test_id AS hi_n_test_id,
    HINTAR.testedness AS hint_testedness,
    HINTAR.n_test_id AS hint_n_test_id,
    CASE
      WHEN HIAR.avg_logfold <= -3 THEN "LOW"
      WHEN HIAR.avg_logfold > -3 AND HIAR.avg_logfold < -1 THEN "LIKE(ISH)"
      WHEN HIAR.avg_logfold >= -1 THEN "LIKE"
      ELSE NULL
    END AS hi_aggregate_label,
    CASE
      WHEN HINTAR.avg_logfold <= -3 THEN "LOW"
      WHEN HINTAR.avg_logfold > -3 AND HINTAR.avg_logfold < -1 THEN "LIKE(ISH)"
      WHEN HINTAR.avg_logfold >= -1 THEN "LIKE"
      ELSE NULL
    END AS hint_aggregate_label
  FROM
    LINK_AA_SEQUENCE AS LAAS
  LEFT JOIN
    LINK_HI_ANTIGENIC_RECORD AS HIAR ON LAAS.variant_hash = HIAR.variant_hash
  LEFT JOIN
    LINK_HINT_ANTIGENIC_RECORD AS HINTAR ON LAAS.variant_hash = HINTAR.variant_hash
  WHERE
    (
      NOT HIAR.testedness IS NULL
      OR NOT HINTAR.testedness IS NULL
    )
),

-- Merge all variant metadata into a consolidated dataset
-- Exclude all VARIANT_HASH entries with archived unrelaxed PDB models and/or relaxed PDB model metrics/parameters
CONSOLIDATED_METADATA AS (
  SELECT
    DISTINCT
    LSM.ctype,
    LSM.subtype,
    LSM.genetic_group,
    LSM.clade,
    LSM.protein,
    LSM.variant_hash,
    LSM.rep_nt_id,
    LSM.cumulative_surveillance,
    LSM.windowed_surveillance,
    LSM.internal_surveillance,
    LSM.external_surveillance,
    LSM.last_seen,
    IF(
      (
        NOT BAT.hi_testedness IS NULL
        OR NOT BAT.hint_testedness IS NULL
      ),
      TRUE,
      FALSE
    ) AS antigenic_record,
    BAT.hi_testedness,
    BAT.hi_n_test_id,
    BAT.hint_testedness,
    BAT.hint_n_test_id,
    BAT.hi_aggregate_label,
    BAT.hint_aggregate_label,
    LSM.cds_seq_length,
    LSM.aa_seq_length,
    LSM.shift_indel,
    LSM.has_insertion,
    DSSV.signature_changes,
    LSM.aa_seq,
    LSM.variant_substitutions,
    LSM.variant_base_substitutions
  FROM
    LINK_SURVEILLANCE_METRICS AS LSM
  LEFT JOIN
    DETECT_SIGNATURE_SEQUENCE_VARIANTS AS DSSV ON LSM.variant_hash = DSSV.variant_hash
  LEFT JOIN
    BIN_ANTIGENIC_TRENDS AS BAT ON LSM.variant_hash = BAT.variant_hash
    )
),

-- Prioritize influenza protein variants (i.e., variant_hash) via DENSE_RANK() computation
-- Implemented "bins" in priortized schema:
--    (1) No antigenic testing results available with a combined number of CDC samples and global specimens >= 50
--    (2) Antigenic results available and a combined number of CDC samples and global specimens >= 50 detected
--    (3) Minimal CDC samples available (N < 3) yet notable global detection (N >= 15)
--    (4) CDC samples available (N >= 3) in addition to notable global detection (N >= 25)
--    (5) Remainder ...
-- Prioritization metrics (within above "bins") for relevant protein variants:
--    (1) Prioritize specimens with antigenic results based on LIKE/LOW categorizations (based on seasonal "reference")
--        (1a) LOW [H1/B] and LOW/LOW [H3 - both HI/HINT protocols]
--        (1b) LOW/LIKE(ISH) [H3 - combination across HI/HINT protocols]
--        (1c) LOW [H3 - either HI/HINT protocol]
--        (1d) LIKE(ISH) [H1/B] and LIKE(ISH)/LIKE(ISH) [H3 - both HI/HINT protocols] and N Test IDs >= 3 (either HI/HINT for H3s)
--        (1e) LIKE(ISH) [H3 - in either HI or HINT protocols] AND N Test IDs >= 3 (either HI/HINT for H3s)
--        (1f) NULL entries for all subtype-specific antigenic protocols
--        (1g) LIKE(ISH) [H1/B] and LIKE(ISH)/LIKE(ISH) [H3 - both HI/HINT protocols] and N Test IDs < 3 (either HI/HINT for H3s)
--        (1h) LIKE(ISH) [H3 - in either HI or HINT protocols] AND N Test IDs < 3 (either HI/HINT for H3s)
--        (1i) Remainder ... (includes all "LIKE" designations)
--    (2) Number of all available CDC and global specimens combined (DESC)
--    (3) Number of "signature" substitutions (DESC)
--    (4) Number of available CDC specimens (DESC)
--    (5) Number of available global specimens (without CDC overlap) (DESC)
--    (6) Clade annotation (DESC) // "newer" clade labels are more advanced alphabetically
--    (7) Variant hash identifier (ASC)
MODEL_PRIORITIZATION AS (
  SELECT
    CMETA.*,
    DENSE_RANK() OVER(
      PARTITION BY
        CMETA.subtype
      ORDER BY
        CASE
          WHEN NOT CMETA.antigenic_record
            AND CMETA.cumulative_surveillance >= 50
            THEN 1
          WHEN CMETA.antigenic_record
            AND CMETA.cumulative_surveillance >= 50
            THEN 2
        END ASC,
        CASE
          WHEN CMETA.internal_surveillance < 3
            AND CMETA.cumulative_surveillance >= 15
            THEN 1
          WHEN CMETA.internal_surveillance >= 3
            AND CMETA.cumulative_surveillance >= 25
            THEN 2
          ELSE 3
        END ASC,
        CASE
          WHEN (
              CMETA.subtype NOT IN ('H3')
              AND CMETA.hi_aggregate_label = 'LOW'
            )
            OR (
              CMETA.subtype IN ('H3')
              AND CMETA.hi_aggregate_label = 'LOW'
              AND CMETA.hint_aggregate_label = 'LOW'
            ) THEN 1
          WHEN CMETA.subtype IN ('H3')
            AND (
              (
                CMETA.hi_aggregate_label = 'LOW'
                AND CMETA.hint_aggregate_label = 'LIKE(ISH)'
              )
              OR (
                CMETA.hi_aggregate_label = 'LIKE(ISH)'
                AND CMETA.hint_aggregate_label = 'LOW'
              )
            ) THEN 2
          WHEN CMETA.subtype IN ('H3')
            AND (
              CMETA.hi_aggregate_label = 'LOW'
              OR CMETA.hint_aggregate_label = 'LOW'
            ) THEN 3
          WHEN (
              CMETA.subtype NOT IN ('H3')
              AND CMETA.hi_aggregate_label = 'LIKE(ISH)'
              AND CMETA.hi_n_test_id >= 3
            )
            OR (
              CMETA.subtype IN ('H3')
              AND CMETA.hi_aggregate_label = 'LIKE(ISH)'
              AND CMETA.hint_aggregate_label = 'LIKE(ISH)'
              AND (
                CMETA.hi_n_test_id >= 3
                OR CMETA.hint_n_test_id >= 3
              )
            ) THEN 4
          WHEN CMETA.subtype IN ('H3')
            AND (
              CMETA.hi_aggregate_label = 'LIKE(ISH)'
              OR CMETA.hint_aggregate_label = 'LIKE(ISH)'
            ) 
            AND (
                CMETA.hi_n_test_id >= 3
                OR CMETA.hint_n_test_id >= 3
            ) THEN 5
          WHEN (
              CMETA.subtype NOT IN ('H3')
              AND CMETA.hi_aggregate_label IS NULL
            )
            OR (
              CMETA.subtype IN ('H3')
              AND CMETA.hi_aggregate_label IS NULL
              AND CMETA.hint_aggregate_label IS NULL
            ) THEN 6
          WHEN (
              CMETA.subtype NOT IN ('H3')
              AND CMETA.hi_aggregate_label = 'LIKE(ISH)'
              AND CMETA.hi_n_test_id < 3
            )
            OR (
              CMETA.subtype IN ('H3')
              AND CMETA.hi_aggregate_label = 'LIKE(ISH)'
              AND CMETA.hint_aggregate_label = 'LIKE(ISH)'
              AND (
                CMETA.hi_n_test_id < 3
                AND CMETA.hint_n_test_id < 3
              )
            ) THEN 7
          WHEN CMETA.subtype IN ('H3')
            AND (
              CMETA.hi_aggregate_label = 'LIKE(ISH)'
              OR CMETA.hint_aggregate_label = 'LIKE(ISH)'
            ) 
            AND (
                CMETA.hi_n_test_id < 3
                OR CMETA.hint_n_test_id < 3
            ) THEN 8
          ELSE 9
        END ASC,
        CMETA.cumulative_surveillance DESC,
        CMETA.signature_changes DESC,
        CMETA.internal_surveillance DESC,
        CMETA.external_surveillance DESC,
        CMETA.clade DESC,
        CMETA.variant_hash
    ) AS model_priority
  FROM CONSOLIDATED_METADATA AS CMETA
  WHERE
    CMETA.windowed_surveillance > 0
    AND (
      CMETA.cumulative_surveillance >= 3
      OR CMETA.antigenic_record
    )
)


------------------------
--- SELECT STATEMENT ---
------------------------

INSERT
  OVERWRITE TABLE protein_modeling.model_priority
SELECT
  DISTINCT
    MP.*
FROM
  MODEL_PRIORITIZATION AS MP;egg_passage
