test_annotateWithGeneParts = function(){
    # bed with two features:
    # 10 000 - 20 000 (exons: 12 000 - 13 000, 14 000 - 15 000, 16 000 - 17 000)
    bed12 = system.file("unitTests/bed12.bed", package = "genomation")
    features = suppressMessages(readTranscriptFeatures(bed12))

    # 4 features with respect to above bed12:
    # intergenic
    # promoter
    # intronic
    # exonic
    target = GRanges(
        rep('ch1', 4),
        IRanges(c(8001, 10001, 13001, 14001), c(8500, 10500, 13500, 14500))
    )
    result = annotateWithGeneParts(target, features)

    members_expected = matrix(
        c(0, 0, 0,
          1, 0, 0,
          0, 0, 1,
          0, 1, 0),
        ncol = 3, nrow = 4, byrow = TRUE,
        dimnames = list(NULL, c("prom", "exon", "intron"))
    )
    checkIdentical(members_expected, result@members)


    # distance is calculated (for both positive strand):
    # start/end - TSS - 1
    # if > 0, then +2 (last base of codon)
    # smaller of abs(dist_start), abs(dist_end)
    dist_to_tss_expected = data.frame(
        target.row = c(1L, 2L, 3L, 4L),
        dist.to.feature = c(-1501, 0, 3002, 4002),
        feature.name = rep("f1", 4),
        feature.strand = factor("+", levels = c("+", "-", "*")),
        row.names = make.unique(rep("1", 4))
        )

    checkIdentical(dist_to_tss_expected, result@dist.to.TSS)

    df = as.data.frame(result)
    df_expected = data.frame(
        "prom" = c(0, 1, 0, 0),
        "exon" = c(0, 0, 0, 1),
        "intron" = c(0, 0, 1, 0),
        "dist.to.feature" = c(-1501, 0, 3002, 4002),
        "feature.name" = rep("f1", 4),
        row.names = 1:4
        )
    checkIdentical(df_expected, df)
}
