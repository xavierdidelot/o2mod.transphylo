## ---- echo = FALSE-------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>", 
  fig.width=8, 
  fig.height=5, 
  fig.path="figs-introduction/"
)

## ------------------------------------------------------------------------
library(ape)
library(outbreaker2)
library(TransPhylo)
library(o2mod.transphylo)
set.seed(0)

## ------------------------------------------------------------------------
nsam <- 5
phy <- rtree(nsam, br=rexp)
phy$tip.label <- 1:nsam
phy$edge.length <- round(phy$edge.length * 365)
plot(phy)
axisPhylo()

## ------------------------------------------------------------------------
library(distcrete)
dates <- dist.nodes(phy)[nsam+1,1:nsam]
si <- distcrete("gamma", shape = 2, scale = 0.5 * 365, interval = 1L, w = 0)
w <- si$d(1:(2*max(dates)))
data <- outbreaker_data(dates = dates, w_dens = w, dna = as.DNAbin(matrix('A',nsam,1)))
data$ptree <- ptreeFromPhylo(phy,max(dates))

## ------------------------------------------------------------------------
res=o2mod.transphylo(data)

## ------------------------------------------------------------------------
plot(res)

## ------------------------------------------------------------------------
plot(res,type='alpha',burnin=0.1*length(res$step))

## ------------------------------------------------------------------------
plot(res, type = 't_inf', burnin = 0.1 * max(res$step))

## ------------------------------------------------------------------------
tinf <- rep(0, nsam);
alpha <- rep(0, nsam)
for (i in 1:nsam) {
  alpha[i] <- tail(res[[8 + i]],1)
  tinf[i] <- tail(res[[8 + nsam+i]],1)
}

alpha[which(is.na(alpha))] <- 0
ttree <- list(ttree = cbind(tinf, data$dates, alpha),
              nam = data$ptree$nam)
ctree <- combine(ttree, data$ptree)
plotCTree(ctree)

