library(outbreaker2)
library(TransPhylo)

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

lik_TransPhylo <- function(data, param) {
  ttree <- list(ttree = cbind(param$t_inf, data$dates, param$alpha),
                nam = data$ptree$nam)
  ttree$ttree[which(is.na(ttree$ttree[,3])),3] <- 0
  txt <- utils::capture.output(ctree <- combine(ttree,data$ptree))
  if (length(txt)>0) return(-Inf)
  return(TransPhylo::probPTreeGivenTTree(ctree, neg = 365 * 0.25))
}

deduceAlpha <- function(ptree,dates,alpha,t_inf) {
  ttree <- list(ttree = cbind(t_inf, dates, alpha),nam = ptree$nam)
  ttree$ttree[which(is.na(ttree$ttree[,3])),3] <- 0
  txt <- utils::capture.output(ctree <- combine(ttree,ptree))
  if (length(txt)>0 && nchar(txt)>9) return(alpha)
  alpha2=TransPhylo::extractTTree(ctree)$ttree[,3]
  return(replace(alpha2,alpha2==0,NA))
}

api <- get_cpp_api()
new_model <- custom_likelihoods(genetic = lik_TransPhylo)
new_move_tinf <- function(param, data, list_custom_ll = new_model) {
  for (i in 1:data$N) {
    current_ll <- api$cpp_ll_all(data,param, i = NULL, list_custom_ll)
    modif <- sample(c(-100:-1,1:100), 1)
    param$t_inf[i] <- param$t_inf[i] + modif
    current_alpha <- param$alpha
    param$alpha <- deduceAlpha(data$ptree,data$dates,param$alpha,param$t_inf)
    new_ll <- api$cpp_ll_all(data, param, i = NULL, list_custom_ll)
    if (log(stats::runif(1)) > (new_ll - current_ll)) {
      param$t_inf[i] <- param$t_inf[i] - modif
      param$alpha <- current_alpha
    }
  }
  return(param)
}
new_moves <- custom_moves(t_inf = new_move_tinf)

outbreaker_transphylo<-function(data,n_iter=1000,sample_every=1) {
  init_tree <- c(NA,rep(1,length(data$dates) - 1))
  init_t_inf <- c(min(data$ptree$ptree[,1]) - 1, data$dates[2:length(data$dates)] - 1)

  conf <- outbreaker2::create_config(init_tree = init_tree, init_t_inf = init_t_inf,
                        init_kappa = 1, init_pi = 1,
                        find_import = FALSE, n_iter = n_iter,
                        sample_every = sample_every,
                        move_mu = FALSE, move_pi = FALSE,
                        move_eps = FALSE, move_lambda = FALSE,
                        move_alpha = FALSE, move_swap_cases = TRUE,
                        move_t_inf = TRUE, move_kappa = FALSE)

  res <- outbreaker2::outbreaker(data = data, config = conf,
                    likelihoods = new_model,
                    moves = new_moves)
}
