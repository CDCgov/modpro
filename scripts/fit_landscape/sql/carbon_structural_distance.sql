SELECT
    RSD.protein,
    RSD.subtype,
    RSD.reference_context,
    RSD.variant_hash AS identifier,
    RSD.atomic_context,
    RSD.atomic_distance,
    RSD.site_confidence,
    RSD.range_limit,
    RSD.site_range,
    RSD.sd

FROM protein_modeling.relaxed_structural_distance AS RSD

WHERE
    RSD.subtype IN ('B','H1','H3')
    AND RSD.atomic_context = 'carbon'