SELECT
    DISTINCT
    ASREF.subtype,
    ASREF.strain_name,
    ASREF.passage,
    ASREF.lot

FROM
    ref_sys.archive_serum_ref AS ASREF

WHERE
    -- Retains ONLY cell inoculum-derived ferret antisera
    NOT ASREF.egg_vaccine

    -- Dynamic filter to identify all "reportable" antisera indicated for the past [3] WHO surveillance seasons
    AND IF(
        MONTH(NOW()) <= 4,
        ASREF.season >= (YEAR(NOW()) - 1),
        IF(
            MONTH(NOW()) >= 10,
            CONCAT(
                CAST(ASREF.season AS STRING),
                ASREF.hemisphere
            ) IN (
                CONCAT(CAST(YEAR(NOW()) AS STRING), "NH"),
                CONCAT(CAST(YEAR(NOW()) AS STRING), "SH"),
                CONCAT(CAST((YEAR(NOW()) - 1) AS STRING), "NH")
            ),
            ASREF.season >= (YEAR(NOW()) - 1)
        )
    )