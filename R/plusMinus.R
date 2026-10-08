#' @title Positive and negative associations
#'
#' @description Calculate positive and negative associations based on Schoener similarity
#'
#' @details Faster implementation of the funciton \code{pos.neg.abun.r2d} from ((citation))
#' that also only returns information relevant to testing the association of abundance and
#' different types of connectivity (plus and minus).
#'
#' @param x site by species matrix (sites as rows, species as columns, abundances in cells)
#' @param alpha the significance level
#' @param B number of replicates for the null model permutations
#'
#' @return an \code{ncol(x)} by \code{ncol(x)} matrix of species-species similarities
#'
#' @author Andy Rominger <ajrominger@@gmail.com>
#'
#' @export

plusMinus <- function(x, alpha = 0.05, B = 999) {
    # metacommunity "abundnace"
    metaX <- colSums(x)

    # observed Schoener distance
    sd <- schoener(x)
    sd <- sd[lower.tri(sd)]

    # null Schoener distances
    nulls <- simulate(nullmodel(x, 'r2dtable'), nsim = B)
    nullsd <- lapply(1:dim(nulls)[3], function(i) {
        d <- schoener(nulls[, , i])
        return(d[lower.tri(d)])
    })

    # combine null with observed
    nullsd <- cbind(sd, do.call(cbind, nullsd))

    # probabilities that the observed species-species distances are >= or <= the null
    ppos <- rowMeans(nullsd >= nullsd[, 1])
    pneg <- rowMeans(nullsd <= nullsd[, 1])

    # significantly positive and negative edges
    epos <- .signi(ppos, alpha)
    eneg <- .signi(pneg, alpha)

    # number of significantly positive and negative edges
    npos <- sum(epos)
    nneg <- sum(eneg)

    # edge list
    elist <- cbind(t(combn(1:ncol(x), 2)), epos, eneg)

    # correlation between centrality and abundance, and means
    corMeanPos <- .abundCenCorMean(elist[elist[, 3] == 1, 1:2, drop = FALSE], metaX)
    corMeanNeg <- .abundCenCorMean(elist[elist[, 4] == 1, 1:2, drop = FALSE], metaX)

    # number of species in each network
    vpos <- length(unique(as.vector(elist[elist[, 3] == 1, 1:2])))
    vneg <- length(unique(as.vector(elist[elist[, 4] == 1, 1:2])))

    return(list(all = c(v = ncol(x)),
                pos = c(n = npos, corMeanPos),
                neg = c(n = nneg, corMeanNeg)))
}


# function to extract probabilities <= alpha
.signi <- function(x, alpha) {
    return(as.integer(x <= alpha))
}


# function to calculate correlation between abundance and centrality
# and mean abundances, unweighted and weighted by centrality
.abundCenCorMean <- function(elist, abund) {
    if(nrow(elist) < 3) {
        return(list(v = NA, rho = c(rho = NA_real_), p = NA, m = NA, wm = NA))
    } else {
        cen <- table(c(elist[, 1], elist[, 2]))
        x <- abund[as.integer(names(cen))]

        # correlation output
        # o <- cor.test(x, cen, method = 'spearman')
        o <- cor(x, cen, method = 'spearman')

        # relative abundance
        relx <- x / sum(abund)

        # return: corr coeff, corr p-val, mean abund, mean abund weighted by centrality
        return(list(v = length(cen),
                    # rho = o$estimate, p = o$p.value,
                    rho = c(rho = o), p = NA,
                    m = mean(relx), wm = weighted.mean(relx, cen)))
    }
}
