.lr.fine <- function(data, cand, grid.n=32, p.offset=0.025, max.iter=50, max.C=9,
                     max.sC=20, max.len.per.probe=1e6
                     ) {
    rC <- .grid.rC(data$seg, data$lr, data$sd, data$len, data$nC, max.C)
    ## optimize p
    rets <- rbindlist(mclapply(seq_len(nrow(cand)), function(i) {
        p0 <- cand[i,p0]
        D0 <- cand[i,D0]
        D1 <- cand[i,D1]
        S1 <- cand[i,S1]
        ## setup variables
        iter <- 0
        pi <- p0
        Di <- D1
        Si <- S1
        ps <- c()
        Ds <- c()
        ## iterate till convergence
        status <- "running"
        while(TRUE) {
            iter <- iter + 1
            ## optimize p
            p.lo <- max(pi-p.offset, 0.01)
            p.hi <- min(pi+p.offset, 0.99)
            pD <- cbind(p=seq(p.lo, p.hi, length.out=grid.n), D=Di)
            optpD <- rbindlist(lapply(seq_len(nrow(pD)), function(i) {
                .llik.rC.p.D.best(rC, pD[i,1], pD[i,2], max.sC, max.len.per.probe)
            }))
            pj <- optpD[,p0]
            Dj <- optpD[,D1]
            Sj <- optpD[,S1]
            if (
                ((abs(pj-pi)<5e-3) && (abs(Dj-Di)<5e-2)) ||
                ((pj %in% ps) && (Dj %in% Ds))) {
                status <- "converged"
            }
            else if (pj<=0.01+1e-9 || pj>=0.99-1e-9 || Dj<1+1e-9 || Dj>6-1e-9) {
                status <- "diverged"
            }
            else if (iter==max.iter) {
                status <- "maxiter"
            }
            else if (Sj-S1>0.02) {
                status <- "aborted"
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
    return(rets)
}
