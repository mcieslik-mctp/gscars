lrCand <- function(grid, opts) {
    x <- as.matrix(dcast.data.table(grid, D0~p0, value.var="L1")[,-1])
    r <- raster(x)
    ## nearest odd integer >= to 3
    rres <- max(2*floor((nrow(r)*opts$res)/2)+1, 3)
    cres <- max(2*floor((ncol(r)*opts$res)/2)+1, 3)
    wind <- matrix(1, nrow=rres, ncol=cres)
    localmax <- focal(r, fun = max.na.rm, w = wind, pad=TRUE, padValue=NA)
    cand <- grid[Which(localmax==r, cells=TRUE)]
    return(cand)
}
