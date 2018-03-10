.addCoverage <- function(gt, t.bam, n.bam) {
    ## normalize, smooth, and gc-correct
    t.cov <- .runMosdepthTile(t.bam, gt)
    n.cov <- .runMosdepthTile(n.bam, gt)
    mx.cov <- cbind(t.cov, n.cov)
    mx.cov <- t(t(mx.cov)/(colSums(mx.cov)/1e6))
    mcols(gt) <- cbind(mcols(gt), cbind(mx.cov, lr.raw=log2(mx.cov[,1]/mx.cov[,2])))

    ## GC correction
    ## if (gc.correct) {
    ##     weight <- ifelse(gt$blacklist==0 & gt$unmasked>0.9, 1, 0)
    ##     gc.residuals <- limma::loessFit(y=gt$lr, x=gt$gc, weight=weight)$residuals
    ##     lr.offset <- lm(gc.residuals~gt$lr, weights=weight)$coefficients[1]
    ##     gt$lr <- gc.residuals-lr.offset
    ## }
    
    ## remove gross outliers
    gt$lr <- ifelse(is.finite(gt$lr.raw), gt$lr.raw, NA_real_)
    if (any(gt$target)) {
        gt$lr[gt$target] <- .smooth.outliers.gr(gt[gt$target], "lr")
    }
    if (any(!gt$target)) {
        gt$lr[!gt$target] <- .smooth.outliers.gr(gt[!gt$target], "lr")
    }
    return(gt)
}

.annotateTiles <- function(genome.tile, seqi, skip.chr) {
    
    bl1.fn <- system.file("extdata/sv-blacklist-10x-hg38-ucsc.bed", package="gscars")
    bl2.fn <- system.file("extdata/hg38.blacklist.bed.gz", package="gscars")
    cyto.fn <- system.file("extdata/hg38.cytoband.bed", package="gscars")
    strict.fn <- system.file("extdata/1000G-strict-unmask-hg38.bed", package="gscars")

    ## blacklist
    bl.1 <- .robust.import(bl1.fn, seqi)
    bl.2 <- .robust.import(bl2.fn, seqi)
    bl <- reduce(c(granges(bl.1), granges(bl.2)))
    tmp <- findOverlaps(genome.tile, bl)
    tmp <- data.table(
        tile=queryHits(tmp),
        blacklist=width(pintersect(genome.tile[queryHits(tmp)], bl[subjectHits(tmp)]))
    )
    setkey(tmp, tile)
    tmp <- tmp[J(seq_along(genome.tile))]
    tmp[is.na(blacklist), blacklist:=0]
    tmp <- tmp[,.(blacklist=sum(blacklist)),by=tile]
    genome.tile$blacklist <- tmp$blacklist / width(genome.tile)

    ## masking
    strict <- .robust.import(strict.fn, seqi)
    tmp <- findOverlaps(genome.tile, strict)
    tmp <- data.table(
        tile=queryHits(tmp),
        unmasked=width(pintersect(genome.tile[queryHits(tmp)], strict[subjectHits(tmp)]))
    )
    setkey(tmp, tile)
    tmp <- tmp[J(seq_along(genome.tile))]
    tmp[is.na(unmasked), unmasked:=0]
    tmp <- tmp[,.(unmasked=sum(unmasked)),by=tile]
    genome.tile$unmasked <- tmp$unmasked / width(genome.tile)
    
    ## GC content
    tmp <- getSeq(BSgenome.Hsapiens.UCSC.hg38, genome.tile)
    genome.tile$gc <- letterFrequency(tmp, "GC", as.prob=TRUE)[,1]

    ## cytobands
    cyto <- .robust.import(cyto.fn, seqi, skip.chr=skip.chr)
    tmp <- findOverlaps(genome.tile, cyto, select="first")
    genome.tile$cytoband <- cyto[tmp]$name
    genome.tile$arm <- paste0(seqnames(genome.tile), str_sub(genome.tile$cytoband, 1, 1))

    ## 
    genome.tile <- sort(genome.tile)
    return(genome.tile)
}

.getGenomeTiles <- function(tile=10000, skip.chr="chrM") {
    ## tile genome
    seqi <- keepStandardChromosomes(Seqinfo(genome="hg38"))
    genome.tile <- tileGenome(dropSeqlevels(seqi, skip.chr), tilewidth=tile, cut.last.tile.in.chrom=TRUE)
    genome.tile$target <- TRUE
    genome.tile <- .annotateTiles(genome.tile, seqi, skip.chr)
    return(genome.tile)
}

.getTargetTiles <- function(tgt.fn, min.gap=50000, shoulder=300, skip.chr="chrM") {
    ## tile targets
    seqi <- keepStandardChromosomes(Seqinfo(genome="hg38"))
    tgt <- .robust.import(tgt.fn, seqi=seqi)
    tgt$name <- NULL
    tgt$score <- NULL
    tgt$target <- TRUE
    gap <- dropSeqlevels(gaps(tgt), skip.chr, pruning.mode="coarse")
    gap <- gap[strand(gap)=="*"]
    gap <- gap[width(gap)>2*shoulder+min.gap]
    gap <- gap-shoulder
    gap$target <- FALSE
    target.tile <- sort(c(tgt, gap))
    target.tile <- .annotateTiles(target.tile, seqi, skip.chr)
    return(target.tile)
}


importGenomeTile <- function(t.bam, n.bam, ...) {
    tile <- .getGenomeTiles(...)
    tile <- .addCoverage(tile, t.bam, n.bam)
    return(tile)
}

importTargetTile <- function(target, t.bam, n.bam, chr.names="ucsc", ...) {
    if (target=="onco1500-v3") {
        target.fn <- system.file(sprintf("extdata/onco1500-v3-targets-hg38-%s.bed", chr.names), package="gscars")
    } else if (target=="agilent-v4") {
        target.fn <- system.file(sprintf("extdata/agilent-v4-targets-hg38-%s.bed", chr.names), package="gscars")
    }
    else {
        stop("Invalid target.")
    }
    tile <- .getTargetTiles(target.fn, ...)
    tile <- .addCoverage(tile, t.bam, n.bam)
    return(tile)
}
