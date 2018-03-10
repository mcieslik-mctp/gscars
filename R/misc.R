.smooth.outliers <- function(data, chr, pos) {
    print(pos)
    obj <- CNA(data, as.integer(factor(chr)), pos,
              data.type = "logratio",
              sampleid = "sample")
    adj <- smooth.CNA(obj)$sample
    
}

.smooth.outliers.gr <- function(gr, data.col) {
    obj <- CNA(mcols(gr)[[data.col]], as.integer(seqnames(gr)), floor((start(gr)+end(gr))/2),
              data.type = "logratio",
              sampleid = "sample")
    adj <- smooth.CNA(obj)$sample
    return(adj)
}

.robust.import <- function(fn, seqi, skip.chr=NULL) {
    tmp <- import(fn)
    valid.seql <- setdiff(intersect(seqlevels(tmp), seqlevels(seqi)), skip.chr)
    tmp <- keepSeqlevels(tmp, valid.seql, pruning.mode="coarse")
    seqlevels(tmp) <- seqlevels(seqi)
    seqinfo(tmp) <- seqi
    return(tmp)
}
