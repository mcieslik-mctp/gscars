plotSeg <- function(sts, sel.chr=NULL) {
    gr <- unlist(sts$tile)
    tmp <- data.table(
        chr=as.character(seqnames(gr)),
        pos=floor((start(gr)+end(gr))/2),
        l2r=mcols(gr)[["lr"]],
        baf=mcols(gr)[["baf"]],
        seg=mcols(gr)[["seg"]]
    )
    if (!is.null(sel.chr)) {
        tmp <- tmp[chr %in% sel.chr]
    }
    plt.l2r <- ggplot(tmp) + facet_grid(chr~.) +
        aes(x=pos, y=l2r, color=factor(as.integer(seg%%3))) +
        geom_point(size=0.5) +
        coord_cartesian(ylim=c(min(-3, min(tmp$l2r, na.rm=TRUE)),
                               max (3, max(tmp$l2r, na.rm=TRUE)))) +
        scale_color_npg(guide=FALSE) +
        theme_pubr()
    
    plt.baf <- ggplot(tmp) + facet_grid(chr~.) +
        aes(x=pos, y=baf, color=factor(as.integer(seg%%3))) +
        geom_point(size=0.5) +
        coord_cartesian(ylim=c(0,0.5))+
        scale_color_npg(guide=FALSE) +
        theme_pubr()
    
    plt <- grid.arrange(plt.l2r, plt.baf, ncol=1)
    return(plt)
}

plotGrid <- function(grid, cand) {
    plt <- ggplot(llik.0) +
        aes(x=p0, y=D0, fill=ifelse(abs(L1)>max(L1), -max(L1), L1)) +
        geom_tile() +
        scale_fill_gradient2(low="blue", mid="white", high="red",
                             name="likelihood"
                             ) +
        geom_point(data=cand, size=2, color="black") +
        theme_pubr(legend="right")
    return(plt)
}

plotStsC <- function(sts, sel.chr=NULL) {
    gr <- unlist(sts$tile)
    tmp <- data.table(
        chr=as.character(seqnames(gr)),
        pos=floor((start(gr)+end(gr))/2),
        l2r=mcols(gr)[["lr"]],
        baf=mcols(gr)[["baf"]],
        C=mcols(gr)[["C"]]
    )
    if (!is.null(sel.chr)) {
        tmp <- tmp[chr %in% sel.chr]
    }
    plt.l2r <- ggplot(tmp) + facet_grid(chr~.) +
        aes(x=pos, y=l2r, color=factor(C)) +
        geom_point(size=0.5) +
        scale_color_npg() +
        theme_pubr()
    
    plt.baf <- ggplot(tmp) + facet_grid(chr~.) +
        aes(x=pos, y=baf, color=factor(C))+
        geom_point(size=0.5) +
        scale_color_npg(guide=FALSE) +
        theme_pubr()
    
    plt <- grid.arrange(plt.l2r, plt.baf, ncol=1)
    return(plt)
}

