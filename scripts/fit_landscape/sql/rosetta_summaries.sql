SELECT
    DISTINCT RRE.total_score,
    RRE.disulfide_potential,
    RRE.attraction_potential,
    RRE.electrostatic_potential,
    RRE.intraresidue_repulsion,
    RRE.weight_aniso_solvation,
    RRE.repulsive_potential,
    RRE.solvation_potential,
    RRE.hbond_bb_sc,
    RRE.hbond_longrange_bb,
    RRE.hbond_sidechain,
    RRE.hbond_shortrange_bb,
    RRE.omega,
    RRE.yhh_planarity,
    RRE.variant_hash AS identifier,
    RRE.`version`
FROM
    protein_modeling.relaxed_rosetta_energy_rmsd AS RRE
    INNER JOIN protein_modeling.relaxed_structural_distance AS RSD ON RRE.variant_hash = RSD.variant_hash
    INNER JOIN (
        SELECT
            RRE.variant_hash,
            MAX(RRE.`version`) AS version_filter
        FROM
            protein_modeling.relaxed_rosetta_energy_rmsd AS RRE
        GROUP BY
            RRE.variant_hash
    ) AS RRE_FILTER ON RRE.variant_hash = RRE_FILTER.variant_hash
    AND RRE.`version` = RRE_FILTER.version_filter
WHERE
    RSD.subtype IN ('B', 'H1', 'H3')
    AND RSD.protein IN ('HA')