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
##' @aliases cds,GenBankAnnot-method
##' @importMethodsFrom GenomicFeatures cds
##' @examples
##' gb = readGenBank(system.file("sample.gbk", package="genbankr"))
##' cds(gb)
##' exons(gb)
##' genes(gb)
##' @return The expected types, \code{GenomicRanges} for most functions,
##' a \code{DNAStrimgSet} for \code{getSeq}, a \code{GenBankAnnot} for
##' \code{annotations}.
##' @export
setMethod("cds", "GenBankAnnot",
          function(x) x@cds)

##' @docType methods
##' @rdname api-methods
##' @aliases exons,GenBankAnnot-method
##' @importMethodsFrom GenomicFeatures exons
##' @export
setMethod("exons", "GenBankAnnot",
          function(x) x@exons)

##' @docType methods
##' @rdname api-methods
##' @aliases genes,GenBankAnnot-method
##' @importMethodsFrom GenomicFeatures genes
##' @export
setMethod("genes", "GenBankAnnot",
          function(x) x@genes)

##' @docType methods
##' @rdname api-methods
##' @aliases transcripts,GenBankAnnot-method
##' @importMethodsFrom GenomicFeatures transcripts
##' @export
setMethod("transcripts", "GenBankAnnot",
          function(x) x@transcripts)

##' @docType methods
##' @rdname api-methods
##' @aliases cds,GenBankFull-method
##' @export
setMethod("cds", "GenBankFull",
          function(x) cds(annotations(x)))

##' @docType methods
##' @rdname api-methods
##' @aliases exons,GenBankFull-method
##' @export
setMethod("exons", "GenBankFull",
          function(x) exons(annotations(x)))

##' @docType methods
##' @rdname api-methods
##' @aliases genes,GenBankFull-method
##' @export
setMethod("genes", "GenBankFull",
          function(x) genes(annotations(x)))

##' @docType methods
##' @rdname api-methods
##' @aliases transcripts,GenBankFull-method
##' @export
setMethod("transcripts", "GenBankFull",
          function(x) transcripts(annotations(x)))


##' @docType methods
##' @rdname api-methods
##' @export
setGeneric("annotations", function(x) standardGeneric("annotations"))
##' @docType methods
##' @rdname api-methods
##' @aliases annotations,GenBankFull-method
##' @export
setMethod("annotations", "GenBankFull",
          function(x) x@annotations)



##' @docType methods
##' @rdname api-methods
##' @aliases getSeq,GenBankFull-method
##' @importMethodsFrom BSgenome getSeq
##' @export
setMethod("getSeq", "GenBankFull",
          function(x, ...) x@sequence)

##' @docType methods
##' @rdname api-methods
##' @aliases getSeq,GenBankAnnot-method
##' @export
setMethod("getSeq", "GenBankAnnot",
          function(x, ...) stop("This object only contains annotation, no raw sequence is available"))

##' @docType methods
##' @rdname api-methods
##' @aliases getSeq,GenBankFile-method
##' @export
setMethod("getSeq", "GenBankFile",
          function(x, ...) parseGenBank(file = x@resource, seq.only = TRUE, ...))

##' @docType methods
##' @rdname api-methods
##' @aliases getSeq,GBAccession-method
##' @export
setMethod("getSeq", "GBAccession",
          function(x, ...) {
    txt = .getGBfromNuccore(x)
    parseGenBank(text=txt, seq.only = TRUE, ...)
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
##' @aliases accession,GenBankAnnot
##' @export
setMethod("accession", "GenBankAnnot",
          function(x) x@accession)

##' @docType methods
##' @rdname gbk-api
##' @aliases accession,GenBankFull
##' @export
setMethod("accession", "GenBankFull",
          function(x) accession(annotations(x)))

##' @docType methods
##' @title GenBank-annotation specific api methods
##' @rdname gbk-api
##' @aliases vers
##' @export
setGeneric("vers", function(x, ...) standardGeneric("vers"))

##' @docType methods
##' @rdname gbk-api
##' @aliases vers,GenBankAnnot
##' @export
setMethod("vers", "GenBankAnnot",
          function(x) x@version)

##' @docType methods
##' @rdname gbk-api
##' @aliases vers,GenBankFull
##' @export
setMethod("vers", "GenBankFull",
          function(x) vers(annotations(x)))

##' @docType methods
##' @importMethodsFrom GenomicFeatures cdsBy
##' @aliases cdsBy,GenBankAnnot
##' @param by character. Factor to group the resulting GRanges by.
##' @rdname api-methods
##' @export
setMethod("cdsBy", "GenBankAnnot",
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
##' @aliases exonsBy,GenBankAnnot
##' @export              
setMethod("exonsBy", "GenBankAnnot",
          function(x, by = c("tx", "gene")) {
    by = match.arg(by)
    if(by == "tx")
        split(exons(x), exons(x)$transcript_id)
    else
        split(exons(x), exons(x)$gene_id)
})

##' @docType methods
##' @rdname api-methods
##' @aliases cdsBy,GenBankFull
##' @export
setMethod("cdsBy", "GenBankFull",
          function(x, by = c("tx", "gene")) cdsBy(annotations(x), by = by)) 
##' @docType methods
##' @rdname api-methods
##' @aliases exonsBy,GenBankFull
##' @export
setMethod("exonsBy", "GenBankFull",
          function(x, by = c("tx", "gene")) exonsBy(annotations(x), by = by)) 



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


setMethod(show, "GenBankAnnot",
          function(object) {
    cat("GenBank Annotations\n")
    .genbanksum(object)
})


setMethod(show, "GenBankFull",
          function(object) {
    cat("GenBank Annotations with Raw Sequence\n")
    .genbanksum(annotations(object))
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
                object@accession, length(object@sources), sum(sapply(object@sources, width))))
    cat(sprintf("%d genes\n%d transcripts\n%d exons/cds elements\n%d variations\n%d other features\n\n",
                length(genes(object)), length(unique(cds(object)$transcript_id)),
                length(exons(object)), 
                length(object@variations), length(object@other_features)))
}



##' @importMethodsFrom GenomeInfoDb isCircular
##' @rdname api-methods
##' @aliases isCircular,GenBankAnnot
##' @export
setMethod("isCircular", "GenBankAnnot",
          function(x) grepl("circular", x@locus))


##' @rdname api-methods
##' @aliases isCircular,GenBankFull
##' @export
setMethod("isCircular", "GenBankFull",
          function(x) isCircular(annotations(x)))

##' @importMethodsFrom GenomeInfoDb seqinfo
##' @rdname api-methods
##' @aliases seqinfo,GenBankAnnot
##' @export
setMethod("seqinfo", "GenBankAnnot",
          function(x) seqinfo(genes(x)))

##' @rdname api-methods
##' @aliases seqinfo,GenBankFull
##' @export
setMethod("seqinfo", "GenBankFull",
          function(x) seqinfo(genes(annotations(x))))

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
##' @return A \code{GenBankFull} or \code{GenBankAnnot} object.
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
##' @param x A GenBankAnnot or GenBankFull object
##' @return A GRanges for the intergenic regions, defined as regions not
##' overlapping any genes defined in the annotations on either strand.
##' @examples
##' gb = readGenBank(system.file("sample.gbk", package="genbankr"))
##' intergenic(gb)
##' @aliases intergenic,GenBankAnnot-method
##' @export
setMethod("intergenic", "GenBankAnnot",
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
        excl = is.na(nrright) || is.na(nrleft)
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


##' @rdname intergenic
##' @docType methods
##' @aliases intergenic,GenBankFull-method
##' @return A \code{GRanges} object
##' @export
setMethod("intergenic", "GenBankFull",
          function(x) {
    intergenic(annotations(x))
})


##' @title Retrieve variantion features
##' 
##' @description Extract the annotated variants from a GenBankAnnot or GenBankFull object
##' 
##' @rdname variants
##' @docType methods
##' @param x a GenBankAnnot or GenBankFull object
##' @return A VRanges containing the variations annotated in the source file
##' @examples
##' gb = readGenBank(system.file("sample.gbk", package="genbankr"))
##' variants(gb)
##' @export
setGeneric("variants", function(x) standardGeneric("variants"))

##' @rdname variants
##' @aliases variants,GenBankAnnot
##' @export
setMethod("variants", "GenBankAnnot", function(x) x@variations)

##' @rdname variants
##' @aliases variants,GenBankFull
##' @export
setMethod("variants", "GenBankFull", function(x) annotations(x)@variations)
    

##' @title Retrieve 'other' features
##'
##' @description Retrieve  the other features (not covered by a different accessor)
##' from the set of annotations
##' @rdname otherFeatures
##' @docType methods
##' @param x a GenBankAnnot or GenBankFull object
##' @return A GRanges containing the features which don't fall into another
##' category (ie not gene, exon, transcript, cds, or variant) annotated in the
##' source file
##' @examples
##' gb = readGenBank(system.file("sample.gbk", package="genbankr"))
##' otherFeatures(gb)
##' @export
setGeneric("otherFeatures", function(x) standardGeneric("otherFeatures"))

##' @rdname otherFeatures
##' @aliases otherFeatures,GenBankAnnot
##' @export
setMethod("otherFeatures", "GenBankAnnot", function(x) x@other_features)

##' @rdname otherFeatures
##' @aliases otherFeatures,GenBankFull
##' @export
setMethod("otherFeatures", "GenBankFull", function(x)
    annotations(x)@other_features)
