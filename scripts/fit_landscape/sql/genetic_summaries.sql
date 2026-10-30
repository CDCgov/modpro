SELECT
    DISTINCT CASE
        WHEN ALN.isl_ha_type IN ('Hvic', 'HvicNvic') THEN 'B vic'
        WHEN ALN.isl_ha_type = 'H1' THEN 'H1 swl'
        ELSE ALN.isl_ha_type
    END AS subtype,
    ALN.protein,
    TXC.clade,
    ALN.variant_hash AS identifier,
    udx.hamming_distance(INFO.aa_aln, ALN.aa_aln) AS hd,
    udx.pcd(INFO.aa_aln, ALN.aa_aln) AS pcd,
    udx.mutation_list_gly(INFO.aa_aln, ALN.aa_aln) AS substitution_list,
    ALN.aa_seq
FROM
    dais.alignments AS ALN
    INNER JOIN geni.reference_info AS INFO ON ALN.isl_ha_type = INFO.subtype
    AND INFO.protein IN ('HA')
    LEFT JOIN (
        SELECT
            DISTINCT TXB.variant_hash,
            TXB.clade
        FROM
            (
                SELECT
                    TXA.variant_hash,
                    FIRST_VALUE(TXA.clade) OVER(
                        PARTITION BY TXA.variant_hash
                        ORDER BY
                            TXA.n DESC
                    ) AS clade
                FROM
                    (
                        SELECT
                            ALN.variant_hash,
                            COALESCE(
                                TAXA.subclade_alias,
                                TAXA.genetic_group
                            ) AS clade,
                            COUNT() OVER(
                                PARTITION BY ALN.variant_hash,
                                COALESCE(
                                    TAXA.subclade_alias,
                                    TAXA.genetic_group
                                )
                            ) AS n
                        FROM
                            dais.alignments AS ALN
                            LEFT JOIN sand.annotation_taxa AS TAXA ON ALN.nt_id = TAXA.nt_id
                    ) AS TXA
            ) AS TXB
    ) AS TXC ON ALN.variant_hash = TXC.variant_hash -- Limit antigenic data consolidation based on available "relaxed" PDB structure availability/analytics
    INNER JOIN protein_modeling.relaxed_structural_distance AS RSD ON ALN.variant_hash = RSD.variant_hash
    AND RSD.subtype IN ('B', 'H1', 'H3')
WHERE
    ALN.protein IN ('HA')
    AND ALN.isl_ha_type IN ('H1', 'H3', 'Hvic')
    AND ALN.reference_id IN ('CALI07', 'BRISBANE60', 'HK4801')