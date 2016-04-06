
library(genbankr)
library(RUnit)


checkVars = function() {
    raw = parseGenbank("vars.gbk")
    feats = raw$FEATURES
    checkEquals(length(feats), 5L, "Reading variation test did not return 5 features (4 variants and source)")
    gba = genbankr::make_gbannot(raw)
    vr = gba@variations
    checkEquals(length(vr), 4L, "Variation test didn't result in 4 variants")
}

checkJoinCompl = function() {
    raw = parseGenbank("compjoin.gbk")
    feats = raw$FEATURES
    checkEquals(length(feats), 5L, "Reading complement join test did not return 5 features (2 genes, 2 cdss, and source)")
    gba = genbankr::make_gbannot(raw)
    cdss = cds(gba)
    checkEquals(length(cdss), 6L, "Complement join/complement order test did not detect all 6 cds fragments (3 for comp(join()) and 3 for comp(order)).")
    
}
