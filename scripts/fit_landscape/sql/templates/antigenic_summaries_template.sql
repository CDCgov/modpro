SELECT
    A.subtype,
    A.test_protocol,
    A.sr_variant_hash,
    A.lot,
    A.ag_variant_hash AS identifier,
    -- IMPORTANT! Data aggregations in this UNION sub-query parsed based on DISTINCT 'sr_lot' field identifiers
    'lot' AS aggregation,
    CAST(MIN(A.titer_logfold) AS DOUBLE) AS min_logfold,
    CAST(MAX(A.titer_logfold) AS DOUBLE) AS max_logfold,
    CAST(AVG(A.titer_logfold) AS DOUBLE) AS avg_logfold,
    CAST(STDDEV(A.titer_logfold) AS DOUBLE) AS sd_logfold,
    CAST(COUNT(A.titer_logfold) AS INT) AS logfold_count
FROM
    (
        SELECT
            -- Compute aggregate analytics for "biological replicates" (same 'isolate_id' in different 'test_id')
            ANT.subtype,
            ANT.test_protocol,
            ANT.ag_cdc_id,
            ANT.ag_isolate_id,
            ANT.ag_variant_hash,
            ANT.sr_variant_hash,
            ANT.lot,
            AVG(ANT.titer_logfold) AS titer_logfold
        FROM
            (
                SELECT
                    -- Compute aggregate analytics for "technical replicates" (same 'isolate_id' in same 'test_id')
                    ANT.test_subtype AS subtype,
                    ANT.test_protocol,
                    ANT.test_id,
                    ANT.ag_cdc_id,
                    ANT.ag_isolate_id,
                    ANT.ag_variant_hash,
                    ANT.sr_variant_hash,
                    ANT.sr_lot AS lot,
                    AVG(ANT.titer_logfold) AS titer_logfold
                FROM
                    dais.antigenic_tests AS ANT
                    INNER JOIN ref_sys.archive_serum_ref AS ASREF ON ANT.sr_lot = ASREF.lot
                    AND NOT ASREF.egg_vaccine
                    INNER JOIN disc_denorm.specimen AS SPCM ON ANT.ag_cdc_id = SPCM.cdc_id
                    INNER JOIN disc_denorm.isolate AS ISL ON ANT.ag_isolate_id = ISL.isolate_id
                    LEFT JOIN ref_sys.variant_taxonomy AS VT ON ANT.ag_variant_hash = VT.variant_hash -- Detection of consensus-level variants in genetically-linked neuraminidase (NA) sequence(s)
                    LEFT JOIN (
                        SELECT
                            ALN.isolate_id,
                            CAST(
                                (
                                    LENGTH(ALN.aa_aln) - LENGTH(REPLACE(ALN.aa_aln, 'X', ''))
                                ) AS INT
                            ) AS na_variants
                        FROM
                            dais.alignments AS ALN
                        WHERE
                            ALN.protein = 'NA'
                    ) AS ALN ON ANT.ag_sequenced_isolate_id = ALN.isolate_id -- Detection of minor variant populations (5-20% positional incidence) in linked HA/NA sequence(s)
                    LEFT JOIN (
                        SELECT
                            A.isolate_id,
                            COUNT(DISTINCT A.position_disc) AS minor_variants
                        FROM
                            (
                                SELECT
                                    VARF.*,
                                    RIGHT(SUBSTRING(ALNP.aa_aln, 1, VARF.position_aa), 1) AS residue_ref
                                FROM
                                    dais_supp.variants_filtered AS VARF
                                    INNER JOIN dais.aligned_proteins AS ALNP ON VARF.variant_hash = ALNP.variant_hash
                                    AND VARF.reference_id = ALNP.reference_id
                                WHERE
                                    VARF.allele_type = "Minority"
                                    AND VARF.protein IN ('HA', 'NA')
                                    AND VARF.reference_id != "BRISBANE60"
                            ) AS A
                        WHERE
                            A.residue != A.residue_ref
                        GROUP BY
                            A.isolate_id,
                            A.flu_sequence_id
                    ) AS VAR ON ANT.ag_sequenced_isolate_id = VAR.isolate_id
                WHERE
                    SPCM.reason_for_submission = 'surveillance' -- Removes all project/research-based specimens
                    AND (
                        SPCM.strain_name IS NOT NULL
                        AND SPCM.strain_name != ''
                    ) -- Isolates MUST have an attributed strain name designation in DAIS/FluLabDB
                    -- Antigenic surveillance data have been "reported" for the affiliated specimen-level identifier (cdc_id)
                    AND ANT.ag_cdc_id IN (
                        SELECT
                            DISTINCT ISL.cdc_id
                        FROM
                            disc_denorm.isolate AS ISL
                        WHERE
                            ISL.date_report IS NOT NULL
                    )
                    AND ANT.test_subtype IN (
                        'H1 swl',
                        'H3',
                        'B vic',
                        'B yam'
                    ) -- Excludes all non-seasonal subtypes (e.g., H1soiv, H3soiv, H5, and H7)
                    AND (
                        ANT.test_protocol = 'fra_protocol'
                        OR ANT.ag_back_titer >= 4
                    ) -- Isolates tested by the Hemagluttination Inhibition (HI) assay MUST have a back-titer >= 4
                    AND regexp_like(ANT.ag_entry_error, '^\\s*$') -- Removes isolates with logged testing "errors"
                    AND ANT.titer_logfold IS NOT NULL -- Requires a valid 'titer_logfold' result/entry
                    AND NOT ANT.titer_error -- BOOLEAN field is TRUE (no DISC-logged antigenic result error)
                    AND ANT.titer_reportable -- Logged antigenic result MUST be reportable (antigen AND antiserum-levels)
                    AND ANT.ag_type = 'test' -- Removes repeated testing for all "reference" viruses in the applied reference panel(s)
                    AND NOT ANT.ag_passage LIKE '%E%' -- Remove all egg-passaged isolates 
                    AND ANT.ag_variant_hash IS NOT NULL -- Sequence data MUST be available for the isolate-linked HA protein sequence
                    -- Cladistic information MUST be available via LABEL, SAND, and/or Nextclade
                    AND (
                        VT.variant_genetic_group IS NOT NULL
                        AND VT.variant_genetic_group != ''
                    ) --Excluded all tested specimens w/consensus NA protein and minor HA/NA protein polymorphisms
                    --ANY isolates with sequence polymorphisms from the original specimen are excluded
                    AND ANT.ag_isolate_id NOT IN (
                        SELECT
                            IPV.isolate_id
                        FROM
                            ref_sys.isolate_pair_variants AS IPV
                    )
                    AND (
                        VAR.minor_variants IS NULL
                        OR VAR.minor_variants = 0
                    )
                    AND (
                        ALN.na_variants IS NULL
                        OR ALN.na_variants = 0
                    ) -- Included 'test_id' in GROUP BY statement to consolidate "technical replicates" within the same performed test
                GROUP BY
                    ANT.test_subtype,
                    ANT.test_protocol,
                    ANT.test_id,
                    ANT.ag_cdc_id,
                    ANT.ag_isolate_id,
                    ANT.ag_variant_hash,
                    ANT.sr_variant_hash,
                    ANT.sr_lot
            ) AS ANT -- Included BOTH 'ag_cdc_id' (specimen-level) and 'ag_isolate_id' (isolate-level) idnetifiers in GROUP BY statement
            -- Required to consolidate all "biological replicates" performed across multiple tests (e.g., "reference" antigens)
        GROUP BY
            ANT.subtype,
            ANT.test_protocol,
            ANT.ag_cdc_id,
            ANT.ag_isolate_id,
            ANT.ag_variant_hash,
            ANT.sr_variant_hash,
            ANT.lot -- Included SQL UNION w/ similar query to that above to aggregate data for the HINT protocol (modified data exclusions/parameters)
        UNION
        SELECT
            ANT.subtype,
            ANT.test_protocol,
            ANT.ag_cdc_id,
            ANT.ag_isolate_id,
            ANT.ag_variant_hash,
            ANT.sr_variant_hash,
            ANT.lot,
            AVG(ANT.titer_logfold) AS titer_logfold
        FROM
            (
                SELECT
                    ANT.test_subtype AS subtype,
                    ANT.test_protocol,
                    ANT.test_id,
                    ANT.ag_cdc_id,
                    ANT.ag_isolate_id,
                    ANT.ag_variant_hash,
                    ANT.sr_variant_hash,
                    ANT.sr_lot AS lot,
                    AVG(ANT.titer_logfold) AS titer_logfold
                FROM
                    dais.antigenic_tests AS ANT
                    INNER JOIN ref_sys.archive_serum_ref AS ASREF ON ANT.sr_lot = ASREF.lot
                    AND NOT ASREF.egg_vaccine
                    INNER JOIN disc_denorm.specimen AS SPCM ON ANT.ag_cdc_id = SPCM.cdc_id
                    INNER JOIN disc_denorm.isolate AS ISL ON ANT.ag_isolate_id = ISL.isolate_id
                    LEFT JOIN ref_sys.variant_taxonomy AS VT ON ANT.ag_variant_hash = VT.variant_hash
                    LEFT JOIN (
                        SELECT
                            ALN.isolate_id,
                            CAST(
                                (
                                    LENGTH(ALN.aa_aln) - LENGTH(REPLACE(ALN.aa_aln, 'X', ''))
                                ) AS INT
                            ) AS na_variants
                        FROM
                            dais.alignments AS ALN
                        WHERE
                            ALN.protein = 'NA'
                    ) AS ALN ON ANT.ag_sequenced_isolate_id = ALN.isolate_id
                    LEFT JOIN (
                        SELECT
                            A.isolate_id,
                            COUNT(DISTINCT A.position_disc) AS minor_variants
                        FROM
                            (
                                SELECT
                                    VARF.*,
                                    RIGHT(SUBSTRING(ALNP.aa_aln, 1, VARF.position_aa), 1) AS residue_ref
                                FROM
                                    dais_supp.variants_filtered AS VARF
                                    INNER JOIN dais.aligned_proteins AS ALNP ON VARF.variant_hash = ALNP.variant_hash
                                    AND VARF.reference_id = ALNP.reference_id
                                WHERE
                                    VARF.allele_type = "Minority"
                                    AND VARF.protein IN ('HA', 'NA')
                                    AND VARF.reference_id != "BRISBANE60"
                            ) AS A
                        WHERE
                            A.residue != A.residue_ref
                        GROUP BY
                            A.isolate_id,
                            A.flu_sequence_id
                    ) AS VAR ON ANT.ag_sequenced_isolate_id = VAR.isolate_id
                WHERE
                    SPCM.reason_for_submission = 'surveillance'
                    AND (
                        SPCM.strain_name IS NOT NULL
                        AND SPCM.strain_name != ''
                    )
                    AND ANT.test_subtype IN (
                        'H1 swl',
                        'H3',
                        'B vic',
                        'B yam'
                    )
                    AND ANT.test_protocol = 'hint_protocol'
                    AND ANT.test_hint_calculation_method = 'AFLUT' -- Includes ONLY high-throughput HINT protocol output
                    AND regexp_like(ANT.ag_entry_error, '^\\s*$')
                    AND ANT.titer_logfold IS NOT NULL
                    AND NOT ANT.titer_error
                    AND ANT.ag_type = 'test'
                    AND NOT ANT.ag_passage LIKE '%E%'
                    AND ANT.ag_variant_hash IS NOT NULL
                    AND (
                        VT.variant_genetic_group IS NOT NULL
                        AND VT.variant_genetic_group != ''
                    ) --Excluded all tested specimens w/consensus NA protein and minor HA/NA protein polymorphisms
                    --ANY isolates with sequence polymorphisms from the original specimen are excluded
                    AND ANT.ag_isolate_id NOT IN (
                        SELECT
                            IPV.isolate_id
                        FROM
                            ref_sys.isolate_pair_variants AS IPV
                    )
                    AND (
                        VAR.minor_variants IS NULL
                        OR VAR.minor_variants = 0
                    )
                    AND (
                        ALN.na_variants IS NULL
                        OR ALN.na_variants = 0
                    )
                GROUP BY
                    ANT.test_subtype,
                    ANT.test_protocol,
                    ANT.test_id,
                    ANT.ag_cdc_id,
                    ANT.ag_isolate_id,
                    ANT.ag_variant_hash,
                    ANT.sr_variant_hash,
                    ANT.sr_lot
            ) AS ANT
        GROUP BY
            ANT.subtype,
            ANT.test_protocol,
            ANT.ag_cdc_id,
            ANT.ag_isolate_id,
            ANT.ag_variant_hash,
            ANT.sr_variant_hash,
            ANT.lot
    ) AS A -- Data aggregated/summarized at isolate-level
    -- Example//Isolates tested with different passage histories (e.g., S1 vs. S2) for the same specimen counted distinctly
WHERE
    A.ag_variant_hash IN ('INSERTVARIANTHASHLIST')
GROUP BY
    A.subtype,
    A.test_protocol,
    A.ag_variant_hash,
    A.sr_variant_hash,
    lot
UNION
SELECT
    B.subtype,
    B.test_protocol,
    B.sr_variant_hash,
    -- Consolidate all 'sr_lot' identifiers with the same viral inoculum/paired HA protein sequence
    udx.sort_list_unique(GROUP_CONCAT(B.lot, ', '), ', ') AS lot,
    B.ag_variant_hash AS identifer,
    -- IMPORTANT! Data aggregations in this UNION sub-query parsed based on DISTINCT 'sr_variant_hash' field identifiers
    'ha_protein' AS aggregation,
    CAST(MIN(B.titer_logfold) AS DOUBLE) AS min_logfold,
    CAST(MAX(B.titer_logfold) AS DOUBLE) AS max_logfold,
    CAST(AVG(B.titer_logfold) AS DOUBLE) AS avg_logfold,
    CAST(STDDEV(B.titer_logfold) AS DOUBLE) AS sd_logfold,
    CAST(COUNT(B.titer_logfold) AS INT) AS logfold_count
FROM
    (
        SELECT
            D.subtype,
            D.test_protocol,
            D.sr_variant_hash,
            D.inoculated_species,
            -- Included to separate data trends against either ferret- or human antisera
            udx.sort_list_unique(GROUP_CONCAT(D.lot, ', '), ', ') AS lot,
            D.ag_variant_hash,
            D.ag_cdc_id,
            D.ag_isolate_id,
            AVG(D.titer_logfold) AS titer_logfold
        FROM
            (
                SELECT
                    ANT.subtype,
                    ANT.test_protocol,
                    ANT.sr_variant_hash,
                    ANT.lot,
                    ANT.inoculated_species,
                    ANT.ag_variant_hash,
                    ANT.ag_cdc_id,
                    ANT.ag_isolate_id,
                    AVG(ANT.titer_logfold) AS titer_logfold
                FROM
                    (
                        SELECT
                            ANT.test_subtype AS subtype,
                            ANT.test_protocol,
                            ANT.test_id,
                            ANT.sr_variant_hash,
                            ANT.sr_lot AS lot,
                            IF(
                                (
                                    SERUM.animal_id IS NULL
                                    OR ANT.sr_lot LIKE '%HUMAN%'
                                    OR ANT.sr_lot LIKE '%uman%'
                                ),
                                'non-ferret',
                                'ferret'
                            ) AS inoculated_species,
                            ANT.ag_variant_hash,
                            ANT.ag_cdc_id,
                            ANT.ag_isolate_id,
                            AVG(ANT.titer_logfold) AS titer_logfold
                        FROM
                            dais.antigenic_tests AS ANT
                            INNER JOIN ref_sys.archive_serum_ref AS ASREF ON ANT.sr_lot = ASREF.lot
                            AND NOT ASREF.egg_vaccine
                            LEFT JOIN disc_denorm.hiserum AS SERUM ON ANT.sr_lot = SERUM.lot
                            INNER JOIN disc_denorm.specimen AS SPCM ON ANT.ag_cdc_id = SPCM.cdc_id
                            INNER JOIN disc_denorm.isolate AS ISL ON ANT.ag_isolate_id = ISL.isolate_id
                            LEFT JOIN ref_sys.variant_taxonomy AS VT ON ANT.ag_variant_hash = VT.variant_hash
                            LEFT JOIN (
                                SELECT
                                    ALN.isolate_id,
                                    CAST(
                                        (
                                            LENGTH(ALN.aa_aln) - LENGTH(REPLACE(ALN.aa_aln, 'X', ''))
                                        ) AS INT
                                    ) AS na_variants
                                FROM
                                    dais.alignments AS ALN
                                WHERE
                                    ALN.protein = 'NA'
                            ) AS ALN ON ANT.ag_sequenced_isolate_id = ALN.isolate_id
                            LEFT JOIN (
                                SELECT
                                    A.isolate_id,
                                    COUNT(DISTINCT A.position_disc) AS minor_variants
                                FROM
                                    (
                                        SELECT
                                            VARF.*,
                                            RIGHT(SUBSTRING(ALNP.aa_aln, 1, VARF.position_aa), 1) AS residue_ref
                                        FROM
                                            dais_supp.variants_filtered AS VARF
                                            INNER JOIN dais.aligned_proteins AS ALNP ON VARF.variant_hash = ALNP.variant_hash
                                            AND VARF.reference_id = ALNP.reference_id
                                        WHERE
                                            VARF.allele_type = "Minority"
                                            AND VARF.protein IN ('HA', 'NA')
                                            AND VARF.reference_id != "BRISBANE60"
                                    ) AS A
                                WHERE
                                    A.residue != A.residue_ref
                                GROUP BY
                                    A.isolate_id,
                                    A.flu_sequence_id
                            ) AS VAR ON ANT.ag_sequenced_isolate_id = VAR.isolate_id
                        WHERE
                            SPCM.reason_for_submission = 'surveillance'
                            AND (
                                SPCM.strain_name IS NOT NULL
                                AND SPCM.strain_name != ''
                            )
                            AND ANT.ag_cdc_id IN (
                                SELECT
                                    DISTINCT ISL.cdc_id
                                FROM
                                    disc_denorm.isolate AS ISL
                                WHERE
                                    ISL.date_report IS NOT NULL
                            )
                            AND ANT.test_subtype IN (
                                'H1 swl',
                                'H3',
                                'B vic',
                                'B yam'
                            )
                            AND (
                                ANT.test_protocol = 'fra_protocol'
                                OR ANT.ag_back_titer >= 4
                            )
                            AND regexp_like(ANT.ag_entry_error, '^\\s*$')
                            AND ANT.titer_logfold IS NOT NULL
                            AND NOT ANT.titer_error
                            AND ANT.titer_reportable
                            AND ANT.ag_type = 'test'
                            AND NOT ANT.ag_passage LIKE '%E%'
                            AND ANT.ag_variant_hash IS NOT NULL
                            AND (
                                VT.variant_genetic_group IS NOT NULL
                                AND VT.variant_genetic_group != ''
                            ) --Excluded all tested specimens w/consensus NA protein and minor HA/NA protein polymorphisms
                            --ANY isolates with sequence polymorphisms from the original specimen are excluded
                            AND ANT.ag_isolate_id NOT IN (
                                SELECT
                                    IPV.isolate_id
                                FROM
                                    ref_sys.isolate_pair_variants AS IPV
                            )
                            AND (
                                VAR.minor_variants IS NULL
                                OR VAR.minor_variants = 0
                            )
                            AND (
                                ALN.na_variants IS NULL
                                OR ALN.na_variants = 0
                            )
                        GROUP BY
                            ANT.test_subtype,
                            ANT.test_protocol,
                            ANT.test_id,
                            ANT.sr_variant_hash,
                            ANT.sr_lot,
                            inoculated_species,
                            ANT.ag_variant_hash,
                            ANT.ag_cdc_id,
                            ANT.ag_isolate_id
                    ) AS ANT
                GROUP BY
                    ANT.subtype,
                    ANT.test_protocol,
                    ANT.sr_variant_hash,
                    ANT.lot,
                    inoculated_species,
                    ANT.ag_variant_hash,
                    ANT.ag_cdc_id,
                    ANT.ag_isolate_id
            ) AS D
        GROUP BY
            D.subtype,
            D.test_protocol,
            D.sr_variant_hash,
            D.inoculated_species,
            D.ag_variant_hash,
            D.ag_cdc_id,
            D.ag_isolate_id
        UNION
        SELECT
            D.subtype,
            D.test_protocol,
            D.sr_variant_hash,
            D.inoculated_species,
            udx.sort_list_unique(GROUP_CONCAT(D.lot, ', '), ', ') AS lot,
            D.ag_variant_hash,
            D.ag_cdc_id,
            D.ag_isolate_id,
            AVG(D.titer_logfold) AS titer_logfold
        FROM
            (
                SELECT
                    ANT.subtype,
                    ANT.test_protocol,
                    ANT.sr_variant_hash,
                    ANT.lot,
                    ANT.inoculated_species,
                    ANT.ag_variant_hash,
                    ANT.ag_cdc_id,
                    ANT.ag_isolate_id,
                    AVG(ANT.titer_logfold) AS titer_logfold
                FROM
                    (
                        SELECT
                            ANT.test_subtype AS subtype,
                            ANT.test_protocol,
                            ANT.test_id,
                            ANT.sr_variant_hash,
                            ANT.sr_lot AS lot,
                            IF(
                                (
                                    SERUM.animal_id IS NULL
                                    OR ANT.sr_lot LIKE '%HUMAN%'
                                    OR ANT.sr_lot LIKE '%uman%'
                                ),
                                'non-ferret',
                                'ferret'
                            ) AS inoculated_species,
                            ANT.ag_variant_hash,
                            ANT.ag_cdc_id,
                            ANT.ag_isolate_id,
                            AVG(ANT.titer_logfold) AS titer_logfold
                        FROM
                            dais.antigenic_tests AS ANT
                            INNER JOIN ref_sys.archive_serum_ref AS ASREF ON ANT.sr_lot = ASREF.lot
                            AND NOT ASREF.egg_vaccine
                            LEFT JOIN disc_denorm.hiserum AS SERUM ON ANT.sr_lot = SERUM.lot
                            INNER JOIN disc_denorm.specimen AS SPCM ON ANT.ag_cdc_id = SPCM.cdc_id
                            INNER JOIN disc_denorm.isolate AS ISL ON ANT.ag_isolate_id = ISL.isolate_id
                            LEFT JOIN ref_sys.variant_taxonomy AS VT ON ANT.ag_variant_hash = VT.variant_hash
                            LEFT JOIN (
                                SELECT
                                    ALN.isolate_id,
                                    CAST(
                                        (
                                            LENGTH(ALN.aa_aln) - LENGTH(REPLACE(ALN.aa_aln, 'X', ''))
                                        ) AS INT
                                    ) AS na_variants
                                FROM
                                    dais.alignments AS ALN
                                WHERE
                                    ALN.protein = 'NA'
                            ) AS ALN ON ANT.ag_sequenced_isolate_id = ALN.isolate_id
                            LEFT JOIN (
                                SELECT
                                    A.isolate_id,
                                    COUNT(DISTINCT A.position_disc) AS minor_variants
                                FROM
                                    (
                                        SELECT
                                            VARF.*,
                                            RIGHT(SUBSTRING(ALNP.aa_aln, 1, VARF.position_aa), 1) AS residue_ref
                                        FROM
                                            dais_supp.variants_filtered AS VARF
                                            INNER JOIN dais.aligned_proteins AS ALNP ON VARF.variant_hash = ALNP.variant_hash
                                            AND VARF.reference_id = ALNP.reference_id
                                        WHERE
                                            VARF.allele_type = "Minority"
                                            AND VARF.protein IN ('HA', 'NA')
                                            AND VARF.reference_id != "BRISBANE60"
                                    ) AS A
                                WHERE
                                    A.residue != A.residue_ref
                                GROUP BY
                                    A.isolate_id,
                                    A.flu_sequence_id
                            ) AS VAR ON ANT.ag_sequenced_isolate_id = VAR.isolate_id
                        WHERE
                            SPCM.reason_for_submission = 'surveillance'
                            AND (
                                SPCM.strain_name IS NOT NULL
                                AND SPCM.strain_name != ''
                            )
                            AND ANT.test_subtype IN (
                                'H1 swl',
                                'H3',
                                'B vic',
                                'B yam'
                            )
                            AND ANT.test_protocol = 'hint_protocol'
                            AND ANT.test_hint_calculation_method = 'AFLUT'
                            AND regexp_like(ANT.ag_entry_error, '^\\s*$')
                            AND ANT.titer_logfold IS NOT NULL
                            AND NOT ANT.titer_error
                            AND ANT.ag_type = 'test'
                            AND NOT ANT.ag_passage LIKE '%E%'
                            AND ANT.ag_variant_hash IS NOT NULL
                            AND (
                                VT.variant_genetic_group IS NOT NULL
                                AND VT.variant_genetic_group != ''
                            ) --Excluded all tested specimens w/consensus NA protein and minor HA/NA protein polymorphisms
                            --ANY isolates with sequence polymorphisms from the original specimen are excluded
                            AND ANT.ag_isolate_id NOT IN (
                                SELECT
                                    IPV.isolate_id
                                FROM
                                    ref_sys.isolate_pair_variants AS IPV
                            )
                            AND (
                                VAR.minor_variants IS NULL
                                OR VAR.minor_variants = 0
                            )
                            AND (
                                ALN.na_variants IS NULL
                                OR ALN.na_variants = 0
                            )
                        GROUP BY
                            ANT.test_subtype,
                            ANT.test_protocol,
                            ANT.test_id,
                            ANT.sr_variant_hash,
                            ANT.sr_lot,
                            inoculated_species,
                            ANT.ag_variant_hash,
                            ANT.ag_cdc_id,
                            ANT.ag_isolate_id
                    ) AS ANT
                GROUP BY
                    ANT.subtype,
                    ANT.test_protocol,
                    ANT.sr_variant_hash,
                    ANT.lot,
                    inoculated_species,
                    ANT.ag_variant_hash,
                    ANT.ag_cdc_id,
                    ANT.ag_isolate_id
            ) AS D
        GROUP BY
            D.subtype,
            D.test_protocol,
            D.sr_variant_hash,
            D.inoculated_species,
            D.ag_variant_hash,
            D.ag_cdc_id,
            D.ag_isolate_id
    ) AS B
WHERE
    B.ag_variant_hash IN ('INSERTVARIANTHASHLIST')
GROUP BY
    B.subtype,
    B.test_protocol,
    B.sr_variant_hash,
    B.ag_variant_hash