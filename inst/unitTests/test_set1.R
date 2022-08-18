
library(genbankr)
library(RUnit)


test_Vars = function() {
    ## raw = parseGenBank("vars.gbk")
    raw = parseGenBank(system.file("unitTests/vars.gbk", package="genbankr"))
    feats = raw$FEATURES
    checkEquals(length(feats), 5L, "Reading variation test did not return 5 features (4 variants and source)")
    gba = genbankr::make_gbrecord(raw)
    vr = gba@variations
    checkEquals(length(vr), 4L, "Variation test didn't result in 4 variants")
    checkTrue(is(getSeq(gba), "DNAStringSet") && !isCircular(gba), "readGenBank did not handle detecting sequence type and circularity from the LOCUS correctly")
}

test_JoinCompl = function() {
    ## raw = parseGenBank("compjoin.gbk")
    raw = parseGenBank(system.file("unitTests/compjoin.gbk", package="genbankr"))
    feats = raw$FEATURES
    checkEquals(length(feats), 5L, "Reading complement join test did not return 5 features (2 genes, 2 cdss, and source)")
    gba = genbankr::make_gbrecord(raw)
    cdss = cds(gba)
    checkEquals(length(cdss), 6L, "Complement join/complement order test did not detect all 6 cds fragments (3 for comp(join()) and 3 for comp(order)).")

}

test_RetSeqFalse = function() {
    raw = readGenBank(system.file("unitTests/vars.gbk",package="genbankr"),
                      ret.seq=FALSE)
    checkTrue(is(raw, "GenBankRecord"), "Checking that readGenBank is callable with ret.seq=FALSE")
    checkTrue(is.null(getSeq(raw)), "Checking that sequence is returned as NULL with ret.seq=FALSE")
}

test_MultVal = function() {
    res1 = readGenBank(system.file("unitTests/multval.gbk", package="genbankr"))
    vars = variants(res1)
    checkTrue(is(vars$note, "CharacterList"), "Checking that unexpected duplicate columns result in CharacterList (singular feature)")
    res2 = readGenBank(system.file("unitTests/multval.gbk", package="genbankr"))
    vars2 = variants(res2)
    checkTrue(is(vars2$note, "CharacterList"), "Checking that unexpected duplicate columns result in CharacterList (one of many features)")
    checkTrue(is(sources(res2)$db_xref, "CharacterList"), "Checking that known multivalue field (db_xref) results in CharacterList even when not duplicated in any features")
}

test_getSeq = function() {
    fil = GenBankFile(system.file("unitTests/vars.gbk", package="genbankr"))
    sq = getSeq(fil)
    checkTrue(is(sq, "DNAStringSet") && Biostrings::width(sq) == 11820)
}


test_genpept = function() {
    fil = system.file("unitTests/ABD64816.1.gpt", package="genbankr")
    res1 = readGenBank(fil)
    checkTrue(is(getSeq(res1), "AAStringSet"), "Checking that sequence is returned as AAStringSet when given genpept file")
}

test_partial = function() {
    fil = system.file("unitTests/partial.gbk", package="genbankr")
    res1 = readGenBank(fil)

    res2 = readGenBank(fil, partial=FALSE)
    checkIdentical(res1, res2, "Checking that partial=NA and partial=FALSE result in the same object")

    res3 = readGenBank(fil, partial=TRUE)
    checkTrue(length(genes(res3)) == 1, "Checking that partial=TRUE retains a partial feature")
}

test_txdb = function() {
    thing = readGenBank(system.file("unitTests/compjoin.gbk", package="genbankr"))
    tx = makeTxDbFromGenBank(thing)


}

test_locustagusage = function() {
    thing = readGenBank(system.file("unitTests/locus_tag.gbk", package="genbankr"))
    cdss = cds(thing)
    txs = transcripts(thing)
    gns = genes(thing)
    checkTrue(all(grepl("PA14", gns$gene_id)), msg = "checking that locus_tag is used as gene_id when present on genes GRanges")
    checkTrue(all(cdss$gene_id %in% gns$gene_id), msg = "checking that cds GRanges gets the right gene_id when locus_tag is available")
    checkTrue(all(txs$gene_id %in% gns$gene_id), msg = "checking that transcripts GRanges gets the right gene_id when locus_tag is available")
}


test_regression_18 = function() {
    thing = readGenBank(system.file("unitTests/test_t_inaccession.gbk", package="genbankr"))
    checkTrue(identical(accession(thing), "Salvelinus_fonti"))
}

