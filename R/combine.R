#Function to combine a transmission tree and a phylogeny into a colored phylogeny
combine <- function(ttree,ptree) {
  nam <- ttree$nam
  ttree <- ttree$ttree
  ptree <- ptree$ptree
  ttree[,1] <- ttree[,1]-0.01 #avoid having transmission and coalescence at the same time
  if (any(ttree[,1]>=ttree[,2])) {print('infection after sampling');return(NULL)}
  if (min(ttree[,1])>min(ptree[,1])) {print('root incompatibility');return(NULL)}
  for (i in 1:length(nam)) {
    if (ttree[i,3]>0 && ttree[i,1] <= ttree[ttree[i,3],1]) {
      print('transmission before infection')
      return(NULL)
    }
  }
  tree <- ptree
  n <- ceiling(nrow(tree)/2)
  tree <- rbind(tree,matrix(0,n,3))
  source <- which(ttree[,3]==0)
  tree[nrow(tree),1] <- ttree[source,1]
  tree[nrow(tree),2] <- 2*n-1
  notsource <- setdiff(1:n,source)
  i2 <- 0
  for (i in notsource) {
    i2 <- i2+1
    f <- which(tree[,2]==i|tree[,3]==i)
    fi <- i
    while (tree[f,1]>ttree[i,1]) {
      fi <- f
      f <- which(tree[,2]==f|tree[,3]==f)
    }
    if (tree[f,2]==fi) tree[f,2]=2*n-1+i2 else tree[f,3]=2*n-1+i2
    tree[2*n-1+i2,2]=fi
    tree[2*n-1+i2,1]=ttree[i,1]
  }
  MySort <- sort(tree[(n+1):nrow(tree),1],decreasing=TRUE,index.return = TRUE)
  ind <- MySort$ix
  for (i in (n+1):nrow(tree)) {
    for (j in 2:3) {
      if (tree[i,j]>n) tree[i,j] <- n + which(ind==tree[i,j]-n)
    }
  }
  tree <- tree[c(1:n,n+ind),]
  tree <- cbind(tree,TransPhylo:::.computeHost(tree))#note access to hidden internal function
  ctree=list(ctree = tree, nam = nam)
  ttree2=TransPhylo::extractTTree(ctree)
  if (any(ttree2$ttree[,3]!=ttree[,3])) {print('inc')}
 return(ctree)
}


