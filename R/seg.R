.absMedDiff <- function(x, y) {
    abd <- abs(median(x, na.rm=TRUE)-median(y, na.rm=TRUE))
    return(abd)
}

.penalty <- function(x, sd.penalty) {
    sd.penalty+1/x
}

.jointSegArm <- function(arm, K, min.seg, min.lr.diff, min.baf.diff, sd.penalty) {
    ## make sure we have enough points to segment
    max.K <- sum(!is.na(arm$lr) & !is.na(arm$baf))/10
    K <- min(K, max.K)
    if (K>1) {
        ## initial segmentation
        seg0 <- jointSeg(cbind(arm$lr, arm$baf), K=K)$bestBkp
        ## skip short segments
        bpt0 <- c(0, seg0, length(arm)) + 1
        len0 <- diff(bpt0)
        min0 <- pmin(len0[-1], len0[-length(len0)])
        seg1 <- seg0[min0>=min.seg]
        if (length(seg1)>1) {
            bpt1 <- c(1, seg1 + 1)
            len1 <- diff(c(bpt1, length(arm)+1))
            idx1 <- rep(seq_along(bpt1), len1)
            ## compute breakpoints stats
            lr1 <- split(arm$lr, idx1)
            baf1 <- split(arm$baf, idx1)
            stat1 <- data.table(
                seg1=seg1,
                lr.diff =sapply(2:length(bpt1), function(i) .absMedDiff( lr1[[i]],  lr1[[i-1]])),
                baf.diff=sapply(2:length(bpt1), function(i) .absMedDiff(baf1[[i]], baf1[[i-1]])),
                min.len =sapply(2:length(bpt1), function(i) min(length(lr1[[i]]), length(lr1[[i-1]])))
            )
            stat1[is.na(lr.diff), lr.diff:=0]
            stat1[is.na(baf.diff), baf.diff:=0]
            stat1[,penalty:=.penalty(min.len, sd.penalty)]
            
            ## filtered segmentation
            seg2 <- stat1[(
                 lr.diff >  min.lr.diff*penalty |
                 baf.diff > min.baf.diff*penalty
            ), seg1]
            bpt2 <- c(1, seg2 + 1)
            len2 <- diff(c(bpt2, length(arm)+1))
            idx2 <- rep(seq_along(bpt2), len2)
            arm$seg <- idx2
        } else {
            arm$seg <- 1L
        }
    } else {
        arm$seg <- 1L
    }
    return(arm)
}

.addJointSeg <- function(gt, ...) {
    gt <- sort(unname(unlist(endoapply(split(gt, gt$arm), .jointSegArm, ...))))
    ## provide globally unique ids
    tmp <- paste(gt$arm, gt$seg)
    gt$seg <- as.integer(factor(tmp, levels=unique(tmp)))
    return(gt)
}


.addBaf <- function(gt, snp, shoulder) {
    hits <- findOverlaps(snp, gt, maxgap = shoulder-1)
    hits <- hits[gt[subjectHits(hits)]$target] # prefer assignment to target
    hits <- hits[!duplicated(queryHits(hits))] # if snp close to two targets pick one
    bad=ifelse(snp[queryHits(hits)]$t.AF<0.5,
               round(   snp[queryHits(hits)]$t.AF  * snp[queryHits(hits)]$t.DP),
               round((1-snp[queryHits(hits)]$t.AF) * snp[queryHits(hits)]$t.DP))
    tmp <- data.table(
        idx=subjectHits(hits),
        bad=bad,
        depth=snp[queryHits(hits)]$t.DP
    )
    setkey(tmp, idx)
    tmp <- tmp[J(1:length(gt))]
    tmp <- tmp[,.(baf=sum(bad)/sum(depth)),by=idx]
    gt$baf <- tmp$baf
    gt$baf <- .smooth.outliers.gr(gt, "baf")
    return(gt)
}

.jointSeg <- function(tile, snp, shoulder, K, min.seg, min.lr.diff, min.baf.diff, sd.penalty) {
    tile <- .addJointSeg(tile, K, min.seg, min.lr.diff, min.baf.diff, sd.penalty)
    tile <- split(tile, tile$seg)
    ## create sts
    seg <- unname(sort(unlist(reduce(tile))))
    snp <- split(snp, findOverlaps(snp, seg, select="first"))
    sts <- list(seg=seg, tile=tile, snp=snp)
    return(sts)
}

segmentGenome <- function(tile, var, min.seg=20, K=100, sd.penalty=1) {
    snp <- filterGenomeGermlineHets(var)
    tile <- .addBaf(tile, snp, 0)
    if (is.null(min.sd.diff)) {
        min.sd.diff <- estimateSd(tile$lr)
    }
    if (is.null(min.baf.diff)) {
        min.baf.diff <- estimateSd(tile$baf)
    }
    sts <- .jointSeg(tile, snp, 0, K, min.seg, min.sd.diff, min.baf.diff, sd.penalty)
    return(sts)
}

segmentTarget <- function(tile, var, shoulder=300, min.seg=3, K=30, sd.penalty=1) {
    snp <- filterTargetGermlineHets(var, tile, shoulder)
    tile <- .addBaf(tile, snp, shoulder)
    if (is.null(min.sd.diff)) {
        min.sd.diff <- estimateSd(tile$lr)
    }
    if (is.null(min.baf.diff)) {
        min.baf.diff <- estimateSd(tile$baf)
    }
    sts <- .jointSeg(tile, snp, shoulder, K, min.seg, min.sd.diff, min.baf.diff, sd.penalty)
    return(sts)
}
