setOldClass("NULL")
setClassUnion("DNAStringSetOrNULL", c("DNAStringSet", "NULL"))

##' @title GenBank data objects
##'
##' 
##' @description These objects represent GenBank annotations 
##' 
##' @rdname GenBank-classes
##' @docType methods
##' @examples
##' gb = readGenBank(system.file("sample.gbk", package="genbankr"))
##' gb
##' @exportClass GenBankRecord
setClass("GenBankRecord", slots = list(genes = "GenomicRanges", cds = "GenomicRanges",
                                      exons = "GenomicRanges",
                                      transcripts = "GenomicRanges",
                                      variations = "VRanges",
                                      sources = "GenomicRanges",
                                      other_features = "GenomicRanges",
                                      locus = "character",
                                      definition = "character",
                                      accession = "character",
                                      version = "character",
                                      source = "ANY",
                                      sequence = "DNAStringSetOrNULL"
                                      )
         )



##' @title GenBank File
##'
##' @description A resource class for use within the rtracklayer framework
##'
##' @rdname gbkfile
##' @docType methods
##' @examples
##' fil = GenBankFile(system.file("sample.gbk", package="genbankr"))
##' gb = import(fil)
##' @exportClass GenBankFile
setClass("GenBankFile", contains = "RTLFile")

##' @rdname gbkfile
##' @docType methods
##' @exportClass GBKFile
##' @aliases GBKFile-class
setClass("GBKFile", contains = "GenBankFile")

##' @rdname gbkfile
##' @docType methods
##' @exportClass GBKFile
##' @aliases GBFile-class
setClass("GBFile", contains = "GenBankFile")

##' @title GBAccession ID class
##' 
##' @description A class representing the (versioned) GenBank accession
##'
##' @rdname GBAccession
##' @examples
##' id = GBAccession("U49845.1")
##' \dontrun{gb = readGenBank(id)}
##' @exportClass GBAccession
setClass("GBAccession", contains="character")

##' @rdname GBAccession
##' @param id A versioned GenBank Accession id
##' @return a \code{GBAccession} object.
##' @export
GBAccession = function(id) {
    new("GBAccession", id)
}
