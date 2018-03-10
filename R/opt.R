## data goptimizeeneration
.llik.data <- function(sts, local.sd) {
    tmp <- data.table(
        seg=unlist(sts$tile)$seg,
        lr=unlist(sts$tile)$lr
    )
    tmp[,":="(n.lr=sum(!is.na(lr))), by=seg]
    global.sd <- estimateSd(tmp$lr)
    if (local.sd) {
        tmp[,sd := sd(lr, na.rm=TRUE), by=seg]
        tmp[,sd := ifelse(n.lr>30, sd, global.sd)]
    } else {
        tmp[,sd := global.sd]
    }
    return(tmp)
}

##
.cand.lr <- function(cand, data) {
    seg <- data[,.(lr=mean(lr, na.rm=TRUE),
                   n.lr=seg.n[1]
                   ), by=seg]
    setkey(cand, seg)
    setkey(seg, seg)
    cand <- cand[seg]
    return(cand)
}

## grid generation
.grid.pD <- function(grid.n, p.lo=0.05, p.hi=0.95, D.lo=1, D.hi=6) {
    p <- seq(p.lo, p.hi, length.out=grid.n)
    D <- seq(D.lo, D.hi, length.out=grid.n)
    pD <- as.matrix(expand.grid(p=p, D=D))
    return(pD)
}

.grid.rC <- function(seg, lr, sd) {
    rC <- expand.grid(lr=lr, C=0:7)
    rC$seg <- seg
    rC$sd <- sd
    setDT(rC)
    setkey(rC, C, seg)
    return(rC)
}

.grid.pDC <- function(cand) {
    pDC <- as.matrix(dcast.data.table(cand, p+D~seg, value.var="C"))
    return(pDC)
}

## lr likelihood function
.llik.rC.p.D <- function(rC, p, D) {
    rC[, ## normal
       ":="(
           norm = dnorm(lr, mean=log2((p*C+(1-p)*2)/D), sd=sd, log=TRUE),
           unif = dunif(lr, min=-6, max=6, log=TRUE)*2
       )
       ]
    x <- rC[,.( ## sum by C and seg
        norm = sum(norm, na.rm=TRUE),
        unif = sum(unif, na.rm=TRUE),
        n = .N
    ), by=.(C, seg)]
    return(x)
}

## 2D grid screening
.llik.grid.inner <- function(rC, p, D) {
    ## compute likelihood for each segment and each C
    x <- .llik.rC.p.D(rC, p, D)
    ## for each segment pick C with highest llik
    x <- x[order(-norm),.SD[1],by=seg]
    ## 
    x[,subc:=(unif>norm)]
    ## result
    x <- list(
        L1 = x[(!subc),sum(norm)], ## segment likelihood
        S1 = x[,weighted.mean(subc, n)], 
        D1 = x[,p * weighted.mean(C, n) + (1 - p) * 2],
        p0 = p, D0 = D
    )
    return(x)
}

.llik.grid <- function(data, grid.n=NULL, pD=NULL) {
    rC <- .grid.rC(data$seg, data$lr, data$sd)
    if (is.null(pD) & !is.null(grid.n)) {
        pD <- .grid.pD(grid.n)
    } else {
        stop("specify grid.n or pD")
    }
    pD <- rbindlist(mclapply(seq_len(nrow(pD)), function(i) {
        .llik.grid.inner(rC, pD[i,1], pD[i,2])
    }, mc.cores=detectCores()))
    return(pD)
}

.llik.cand <- function(llik, res=0.1) {
    x <- as.matrix(dcast.data.table(llik, D0~p0, value.var="L1")[,-1])
    r <- raster(x)
    ## nearest odd integer >= to 3
    rres <- max(2*floor((nrow(r)*res)/2)+1, 3)
    cres <- max(2*floor((ncol(r)*res)/2)+1, 3)
    wind <- matrix(1, nrow=rres, ncol=cres)
    maxf <- function(x) max(x, na.rm=TRUE)
    localmax <- focal(r, fun = maxf, w = wind, pad=TRUE, padValue=NA)
    cand <- llik[Which(localmax==r, cells=TRUE)]
    return(cand)
}

.llik.fine <- function(data, cand, grid.n=32, p.offset=0.025, max.iter=25, group.digits=2) {
    rC <- .grid.rC(data$seg, data$lr, data$sd)
    ## optimize p
    rets <- rbindlist(mclapply(seq_len(nrow(cand)), function(i) {
        p0 <- cand[i,p0]
        D0 <- cand[i,D0]
        D1 <- cand[i,D1]
        ## setup variables
        iter <- 0
        pi <- p0
        Di <- D1
        ps <- c()
        Ds <- c()
        ## iterate till convergence
        status <- "running"
        while(TRUE) {
            iter <- iter + 1
            p.lo <- max(pi-p.offset, 0)
            p.hi <- min(pi+p.offset, 1)
            pD <- cbind(p=seq(p.lo, p.hi, length.out=grid.n), D=Di)
            pD <- rbindlist(mclapply(seq_len(nrow(pD)), function(i) {
                .llik.grid.inner(rC, pD[i,1], pD[i,2])
            }))
            optpD <- pD[order(-L1),.SD[1]]
            pj <- optpD[,p0]
            Dj <- optpD[,D1]
            if (
                ((abs(pj-pi)<1e-3) && (abs(Dj-Di)<1e-3)) ||
                ((pj %in% ps) && (Dj %in% Ds))) {
                status <- "converged"
            }
            else if (abs(pj-p0) > 0.2 || abs(Dj-D1) > 1) {
                status <- "diverged"
            }
            else if (iter==max.iter) {
                status <- "maxiter"
            }
            ps <- c(ps, pj)
            Ds <- c(Ds, Dj)
            pi <- pj
            Di <- Dj
            cat(sprintf("candidate:%s iter:%s p:%.4f D:%.4f llik:%.0f status:%s\n",
                        i, iter, pi, Di, optpD[,L1], status))
            if (status != "running") {
                break
            }
        }
        ret <- list(p0=p0, D0=D0, D1=D1, L1=cand[i,L1], S1=cand[i,S1], pi=pi, Di=Di,
                    Li=optpD[,L1], Si=optpD[,S1], iter=iter, status=status)
        return(ret)
    }, mc.cores=detectCores()))
    rets[,group:=.GRP, by=.(round(pi, group.digits), round(Di, group.digits))]
    return(rets)
}

## candidate solutions
.llik.full.inner <- function(rC, pi, Di) {
    ## compute likelihood for each segment and each C
    x <- .llik.rC.p.D(rC, pi, Di)
    ## for each segment pick C with highest llik
    x <- x[order(-norm),.SD[1],by=seg]
    x[,":="(pi=pi, Di=Di)]
    setkey(x, seg)
    return(x)
}

.llik.full <- function(data, pi, Di) {
    rC <- .grid.rC(data$seg, data$lr, data$sd)
    CL <- rbindlist(mclapply(seq_len(length(pi)), function(i) {
        .llik.full.inner(rC, pi[i], Di[i])
    }, mc.cores=detectCores()))
    CL <- cbind(CL, data[,.(lr=mean(lr,na.rm=TRUE), n.lr=n.lr[1]), by=seg])
    CL[,sC := (2^(lr) * Di)/pi - ((2 * (1 - pi))/pi)]
    return(CL)
}
