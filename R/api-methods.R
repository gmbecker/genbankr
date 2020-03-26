##' @import BiocGenerics
##' @import S4Vectors
NULL



##' @title Annotation extraction api
##'
##'
##' @description Accessor functions shared with the larger Bioconductor
##' ecosystem.
##' 
##' @param x The object containing the annotations
##' @param ... unused.
##' @docType methods
##' @rdname api-methods
##' @aliases cds,GenBankRecord-method
##' @importMethodsFrom GenomicFeatures cds
##' @examples
##' gb = readGenBank(system.file("sample.gbk", package="genbankr"))
##' cds(gb)
##' exons(gb)
##' genes(gb)
##' @return The expected types, \code{GenomicRanges} for most functions,
##' a \code{DNAStrimgSet} for \code{getSeq}
##' @export
setMethod("cds", "GenBankRecord",
          function(x) x@cds)

##' @docType methods
##' @rdname api-methods
##' @aliases exons,GenBankRecord-method
##' @importMethodsFrom GenomicFeatures exons
##' @export
setMethod("exons", "GenBankRecord",
          function(x) x@exons)

##' @docType methods
##' @rdname api-methods
##' @aliases genes,GenBankRecord-method
##' @importMethodsFrom GenomicFeatures genes
##' @export
setMethod("genes", "GenBankRecord",
          function(x) x@genes)

##' @docType methods
##' @rdname api-methods
##' @aliases transcripts,GenBankRecord-method
##' @importMethodsFrom GenomicFeatures transcripts
##' @export
setMethod("transcripts", "GenBankRecord",
          function(x) x@transcripts)

##' @importMethodsFrom Biostrings getSeq
##' @docType methods
##' @rdname api-methods
##' @aliases getSeq,GenBankRecord-method
##' @export
setMethod("getSeq", "GenBankRecord",
          function(x, ...) x@sequence)


##' @docType methods
##' @rdname api-methods
##' @aliases getSeq,GenBankFile-method
##' @export
setMethod("getSeq", "GenBankFile",
          function(x, ...) parseGenBank(file = x@resource, ret.anno=FALSE, ...))

##' @docType methods
##' @rdname api-methods
##' @aliases getSeq,GBAccession-method
##' @export
setMethod("getSeq", "GBAccession",
          function(x, ...) {
    txt = .getGBfromNuccore(x)
    parseGenBank(text=txt, ret.anno = FALSE, ...)
})


##' @docType methods
##' @name gbk-specific-api
##' @title genbankr specific api
##'
##' @description Accessor functions specific to genbankr objects.
##' 
##' @rdname gbk-api
##' @param x A genbank annotation object
##' @param ... unused.
##' @aliases accession
##' @examples
##' gb = readGenBank(system.file("sample.gbk", package="genbankr"))
##' accession(gb)
##' vers(gb)
##' @return Character vectors for \code{accession} and \code{vers}
##' @export
setGeneric("accession", function(x, ...) standardGeneric("accession"))

##' @docType methods
##' @rdname gbk-api
##' @aliases accession,GenBankRecord
##' @export
setMethod("accession", "GenBankRecord",
          function(x) x@accession)



##' @docType methods
##' @rdname gbk-api
##' @aliases definition
##' @export
setGeneric("definition", function(x, ...) standardGeneric("definition"))

##' @docType methods
##' @rdname gbk-api
##' @aliases definition,GenBankRecord
##' @export
setMethod("definition", "GenBankRecord",
          function(x) x@definition)


##' @docType methods
##' @rdname gbk-api
##' @aliases locus
##' @export
setGeneric("locus", function(x, ...) standardGeneric("locus"))

##' @docType methods
##' @rdname gbk-api
##' @aliases locus,GenBankRecord
##' @export
setMethod("locus", "GenBankRecord",
          function(x) x@locus)






##' @docType methods
##' @title GenBank-annotation specific api methods
##' @rdname gbk-api
##' @aliases vers
##' @export
setGeneric("vers", function(x, ...) standardGeneric("vers"))

##' @docType methods
##' @rdname gbk-api
##' @aliases vers,GenBankRecord
##' @export
setMethod("vers", "GenBankRecord",
          function(x) x@version)



##' @docType methods
##' @title GenBank-annotation specific api methods
##' @rdname gbk-api
##' @aliases sources
##' @export
setGeneric("sources", function(x, ...) standardGeneric("sources"))

##' @docType methods
##' @rdname gbk-api
##' @aliases sources,GenBankRecord
##' @export
setMethod("sources", "GenBankRecord",
          function(x) x@sources)




##' @docType methods
##' @importMethodsFrom GenomicFeatures cdsBy
##' @aliases cdsBy,GenBankRecord
##' @param by character. Factor to group the resulting GRanges by.
##' @rdname api-methods
##' @export
setMethod("cdsBy", "GenBankRecord",
          function(x, by = c("tx", "gene")) {
    by = match.arg(by)
    if(by == "tx")
        split(cds(x), cds(x)$transcript_id)
    else
        split(cds(x), cds(x)$gene_id)
})

##' @docType methods
##' @rdname api-methods
##' @importMethodsFrom GenomicFeatures exonsBy
##' @aliases exonsBy,GenBankRecord
##' @export              
setMethod("exonsBy", "GenBankRecord",
          function(x, by = c("tx", "gene")) {
    by = match.arg(by)
    if(by == "tx")
        split(exons(x), exons(x)$transcript_id)
    else
        split(exons(x), exons(x)$gene_id)
})


setMethod(showAsCell, "XStringSet", function(object) {
    wds = width(object)
    short = wds < 9
    end1 = rep(3, times = length(wds))
    end1[short] = wds[short]
    st2 = wds - 2
    pst = rep("...", times = length(wds))
    pst2 = rep("", times = length(wds))
    pst2[!short] = substr(object[!short], st2[!short], wds[!short])
    pst[short] = ""
    paste0(substr(object, 1, end1), pst, pst2)
})


setMethod(show, "GenBankRecord",
          function(object) {
    cat("GenBank Annotations\n")
    .genbanksum(object)
})


setMethod(show, "GBAccession",
          function(object) {
    msg = if(length(object) > 6) {
              paste(paste(head(object, 3), collapse=" "),
                    "...", paste(tail(object,3), collapse=" "))
    } else {
        paste(object, collapse = " ")
    }
    
                 
                 
    cat("GenBank Accession Number(s): ", msg, "\n")
})


.genbanksum = function(object) {
        cat(sprintf("%s \nAccession: %s\n%d Sequence(s) with total length length: %d\n", object@definition,
                object@accession, length(object@sources), sum(width(object@sources))))
    cat(sprintf("%d genes\n%d transcripts\n%d exons/cds elements\n%d variations\n%d other features\n\n",
                length(genes(object)), length(unique(cds(object)$transcript_id)),
                length(exons(object)), 
                length(object@variations), length(object@other_features)))
}



##' @importMethodsFrom GenomeInfoDb isCircular
##' @rdname api-methods
##' @aliases isCircular,GenBankRecord
##' @export
setMethod("isCircular", "GenBankRecord",
          function(x) grepl("circular", x@locus))


##' @importMethodsFrom GenomeInfoDb seqinfo
##' @rdname api-methods
##' @aliases seqinfo,GenBankRecord
##' @export
setMethod("seqinfo", "GenBankRecord",
          function(x) seqinfo(genes(x)))

##' @importMethodsFrom rtracklayer import
##' @importFrom rtracklayer resource
##' @title Import genbank file
##' @description Import a genbank file using the rtracklayer API.
##' @docType methods
##' @param con See import docs.
##' @param format See import docs.
##' @param text See import docs.
##' @param ... Arguments passed to \code{readGenBank}
##' @aliases import,GenBankFile
##' @return A \code{GenBankRecord} object.
##' @rdname import
##' @export
setMethod("import", "GenBankFile",
          function(con, format, text, ...) 
              {
                  fil = resource(con)
                  readGenBank(file = fil, ...)
              })
##' @title GenBank file for use with import
##'
##' @description Create a GenBankFile object.
##' 
##' @param fil character. Path to the genbank file
##' @return A \code{GenBankFile} object
##' @rdname gbkfile
##' @export
GenBankFile = function(fil) new("GenBankFile", resource = fil)


##' @importMethodsFrom VariantAnnotation intergenic
##' @title Extract intergenic regions from processed GenBank annotations
##' 
##' @description Extract the intergenic regions from a set of GenBank annotations.
##' 
##' @docType methods
##' @rdname intergenic
##' @name intergenic
##' @param x A GenBankRecord object
##' @return A GRanges for the intergenic regions, defined as regions not
##' overlapping any genes defined in the annotations on either strand.
##' @examples
##' gb = readGenBank(system.file("sample.gbk", package="genbankr"))
##' intergenic(gb)
##' @aliases intergenic,GenBankRecord-method
##' @export
setMethod("intergenic", "GenBankRecord",
          function(x) {
    gns = genes(x)
    strand(gns) = "*"
    spl = split(gns, seqnames(gns))
    res = GRangesList(lapply(spl, function(chrgns) {
        src = x@sources[seqnames(x@sources) == seqnames(chrgns)[1] ]
                       
        ig = gaps(chrgns, start(src), end(src))
        ig = ig[strand(ig) == "*"]
        nrright = precede(ig, chrgns)
        nrleft = follow(ig, chrgns)
        lftlab = chrgns$gene[na.omit(nrleft)]
        rtlab = chrgns$gene[na.omit(nrright)]
        if(is.na(nrleft[1]))
            lftlab = c("SEQ-BEGIN", lftlab)
        if(is.na(nrright[length(nrright)]))
            rtlab = c(rtlab, "SEQ-END")
        
        ig$intergenic_id = paste("intergenic", lftlab, rtlab, sep="_")
        ig
    }))
    unlist(res)
})



##' @title Retrieve variantion features
##' 
##' @description Extract the annotated variants from a GenBankRecord object
##' 
##' @rdname variants
##' @docType methods
##' @param x a GenBankRecord object
##' @return A VRanges containing the variations annotated in the source file
##' @examples
##' gb = readGenBank(system.file("sample.gbk", package="genbankr"))
##' variants(gb)
##' @export
setGeneric("variants", function(x) standardGeneric("variants"))

##' @rdname variants
##' @aliases variants,GenBankRecord
##' @export
setMethod("variants", "GenBankRecord", function(x) x@variations)

##' @title Retrieve 'other' features
##'
##' @description Retrieve  the other features (not covered by a different accessor)
##' from the set of annotations
##' @rdname otherFeatures
##' @docType methods
##' @param x a GenBankRecord object
##' @return A GRanges containing the features which don't fall into another
##' category (ie not gene, exon, transcript, cds, or variant) annotated in the
##' source file
##' @examples
##' gb = readGenBank(system.file("sample.gbk", package="genbankr"))
##' otherFeatures(gb)
##' @export
setGeneric("otherFeatures", function(x) standardGeneric("otherFeatures"))

##' @rdname otherFeatures
##' @aliases otherFeatures,GenBankRecord
##' @export
setMethod("otherFeatures", "GenBankRecord", function(x) x@other_features)


setGeneric("sequence<-", function(x, value) standardGeneric("sequence<-"))
setMethod("sequence<-", "GenBankRecord", function(x, value) {
    x@sequence = value
    x
})
