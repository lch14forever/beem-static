##' @title diagnoseBiomass
##'
##' @param beem.out output of a beem run
##' @param true.biomass measured/true biomass (default: the biomass of last iteration)
##' @param alpha transparency parameter of the biomass lines
##' @description plot the trace of the fitted error of biomass
##' @author Chenhao Li, Niranjan Nagarajan
##' @export
diagnoseBiomass <- function(beem.out, true.biomass=NA, alpha=0.1,...){
    trace.m <- beem.out$trace.m
    nIter <- ncol(trace.m)
    if(any(is.na(true.biomass))){
        true.m <- trace.m[,nIter]
    }else{
        true.m <- true.biomass
    }
    rel.err.m <- (t((trace.m-true.m)/true.m))*100
    col <- rep(rgb(0,0,0,alpha), nrow(trace.m))
    col[beem.out$sample2rm] <- rgb(1,1,0, alpha)
    matplot(rel.err.m, type='l',
            xlab="Iterations",
            ylab="Relative difference (%)",
            main='Biomass trace',
            lty=1, lwd = 3,
            col=col,
            ...
            )
    lines(x=1:ncol(trace.m),y=apply(rel.err.m,1,median), col='red', lwd=5)
}

##' @title diagnoseFit
##'
##' @param beem.out output of a beem run
##' @param dat input data for the beem run
##' @param thre threshold of R2 to show the label
##' @param annotate display the label for the species with R2 > thre
##' @import ggplot2
##' @import ggrepel
##' @description plot the R2 value for each species
##' @author Chenhao Li, Niranjan Nagarajan
##' @export
diagnoseFit <- function(beem.out, dat, thre=0.5, annotate=TRUE){
    dat.tss <- tss(dat)
    r_ss <- rowSums(beem.out$err.p, na.rm=TRUE)
    t_ss <- apply(dat.tss, 1, function(x) sum((x[x!=0]-mean(x[x!=0]))^2))
    r2 <- 1-r_ss/t_ss
    plot.dat <- data.frame(r2=r2, species=rownames(dat.tss))

    p <- ggplot(plot.dat, aes(y=r2, x=species, label=species)) +
        geom_point(size=2) +
        geom_abline(intercept=thre, slope=0, lty=2, size=1.5, col='red') +
        labs(x='Species', y=expression(R^{2})) +
        coord_flip() +
        theme_bw()
    if(annotate){
        p <- p +
            geom_text_repel(data=subset(plot.dat, r2>thre))
    }
    p
}

##' @title showInteraction
##'
##' @param beem.out output of a beem run
##' @param dat input data for the beem run
##' @import igraph
##' @import ggraph
##' @description generate a plot for the interaction network inferred using ggraph
##' @author Chenhao Li, Niranjan Nagarajan
##' @export
showInteraction <- function(beem.out, dat){
    b <- beem2param(beem.out)$b
    diag(b) <- 0
    g <- graph.adjacency(b, mode='directed', weighted='I')
    V(g)$label <- rownames(dat)
    V(g)$RelativeAbundance <- rowMeans(dat)
    E(g)$Type <- ifelse( E(g)$I >0, '+', '-')
    E(g)$Strength <- abs(E(g)$I)
    g.simple <- delete.vertices(g, V(g)[degree(g) == 0])
    ggraph(g.simple, layout = 'fr')+#, circular=TRUE) +
        geom_edge_arc(aes(col=Type, width=Strength),arrow = arrow(length = unit(2, 'mm')),
                      curvature = 0.1, alpha=0.8,
                      end_cap=circle(1.5, 'mm'), start_cap=circle(1.5, 'mm')) +
        geom_node_point(pch=1, aes(size=RelativeAbundance)) +
        geom_node_text(aes(label = label), size=2, repel = TRUE) +
        scale_edge_width(range = c(0.5,1.5), guide=FALSE) +
        theme_void()
}

##' @title pcoa
##'
##' @param countData OTU/species abundance table (each row is one species, each column is one site)
##' @param col a vector of colors for the points
##' @description perform a PCoA analysis using Bray-Curtis distance
##' @importFrom vegan vegdist
##' @import ggplot2
##' @author Chenhao Li, Niranjan Nagarajan
##' @export
pcoa <- function(countData, col='Color'){
    dat.pcoa <- cmdscale(vegan::vegdist(t(countData)), eig=TRUE)
    per.var <- (dat.pcoa$eig/sum(dat.pcoa$eig))[1:2] * 100
    dat <- data.frame(dat.pcoa$points, col=col)
    ggplot(dat, aes(x=X1, y=X2, col=col)) +
        geom_point(size=2) +
        labs(x=paste0('CMD1 (' ,round(per.var[1], 2),'%)'),
             y=paste0('CMD2 (',round(per.var[2], 2),'%)'))
}


##' @title cluster
##'
##' @param countData OTU/species abundance table (each row is one species, each column is one site)
##' @description perform a hierachical clustering analysis using Bray-Curtis distance
##' @importFrom vegan vegdist
##' @author Chenhao Li, Niranjan Nagarajan
##' @export
cluster <- function(countData){
    hc <- hclust(vegan::vegdist(t(countData)))
    hc
}

##' @title auc.b
##'
##' @param b.est The estimated interaction network (This can be a correlation matrix and the edges will be ranked based on the absolute value of the correlation coefficents)
##' @param b.true The true interaction network
##' @param is.association b.est is an association structure
##' @author Chenhao Li, Niranjan Nagarajan
##' @description plot ROC curve with AUC for the interaction graph (interaction signs are ignored)
##' @export
auc.b <- function(b.est, b.true, is.association=FALSE, ...){
    if (!requireNamespace("pROC", quietly = TRUE)) {
        stop("Package \"pROC\" needed for this function to work. Please install it.",
             call. = FALSE)
    }
    diag(b.est) <- NA
    diag(b.true) <- NA
    if(is.association){
        b.true <- interaction2association(b.true)
        est <- b.est[lower.tri(b.est)]
        lab <- b.true[lower.tri(b.true)]
    }else{
        est <- c(b.est)
        lab <- c(b.true)
    }
    est <- abs(est[!is.na(est)])
    lab <- lab[!is.na(lab)]
    lab <- (lab!=0 ) *1
    pROC::plot.roc(lab,est, print.auc = TRUE, ...)
}

##' @title interaction2association
##'
##' @param m an interaction matrix
##' @author Chenhao Li, Niranjan Nagarajan
##' @description Convert interaction matrix to association matrix by taking the average of the symmetric entries
interaction2association <- function(m){
    m.cp <- m
    tmp <- (m[lower.tri(m)] + t(m)[lower.tri(m)])/2
    m.cp[lower.tri(m.cp)] <- tmp
    m.cp[upper.tri(m.cp)] <- tmp
    m.cp
}

