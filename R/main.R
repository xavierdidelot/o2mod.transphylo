#' Run outbreaker2 using the TransPhylo model
#'
#' @param data a list of named items containing input data as returned by
#'     \code{\link[outbreaker2]{outbreaker_data}}, with an additional item
#'     \code{ptree} providing the phylogenetic tree
#'
#' @param config a set of settings as returned by
#'   \code{\link[outbreaker2]{create_config}}
#'
#' @param likelihoods a set of log-likelihood functions as returned by
#'     \code{\link[outbreaker2]{custom_likelihoods}}
#'
#' @param priors a set of log-prior functions as returned by
#'     \code{\link[outbreaker2]{custom_priors}}
#'
#' @param moves a set of movement functions as returned by
#'     \code{\link[outbreaker2]{custom_moves}}
#'
#' @return an \code{outbreaker_chains} object
#'
#' @export
o2mod.transphylo <- function(data,
                             config = NULL,
                             priors = NULL,
                             likelihoods = NULL,
                             moves = NULL) {

  lik_TransPhylo <- function(data, param) {
    ttree <- list(ttree = cbind(param$t_inf, data$dates, param$alpha),
                  nam = data$ptree$nam)
    ttree$ttree[which(is.na(ttree$ttree[, 3])), 3] <- 0
    prob=-Inf
    try({txt <- utils::capture.output(ctree <- combine(ttree, data$ptree))
    if (length(txt) == 0) prob<-TransPhylo::probPTreeGivenTTree(ctree$ctree, neg = 365 * 0.25)
    },silent=T)
    return(prob)
  }

  deduceAlpha <- function(ptree, dates, alpha, t_inf) {
    ttree <- list(ttree = cbind(t_inf, dates, alpha), nam = ptree$nam)
    ttree$ttree[which(is.na(ttree$ttree[, 3])), 3] <- 0
    txt <- utils::capture.output(ctree <- combine(ttree, ptree))
    if (length(txt) > 0 && nchar(txt) > 9) {
      return(alpha)
    }
    alpha2 = TransPhylo::extractTTree(ctree)$ttree[, 3]
    return(replace(alpha2, alpha2 == 0, NA))
  }

  new_move_tinf <- function(param, data, list_custom_ll = def_likelihoods) {
    for (i in 1:data$N) {
      current_ll <- api$cpp_ll_all(data, param, i = NULL, list_custom_ll)
      modif <- sample(c(-100:-1, 1:100), 1)
      param$t_inf[i] <- param$t_inf[i] + modif
      current_alpha <- param$alpha
      param$alpha <- deduceAlpha(data$ptree, data$dates,
                                 param$alpha, param$t_inf)
      new_ll <- api$cpp_ll_all(data, param, i = NULL, list_custom_ll)
      if (log(stats::runif(1)) > (new_ll - current_ll)) {
        param$t_inf[i] <- param$t_inf[i] - modif
        param$alpha <- current_alpha
      }
    }
    return(param)
  }

  init_tree <- c(NA, rep(1, length(data$dates) - 1))
  init_t_inf <- c(min(data$ptree$ptree[, 1]) - 1,
                  data$dates[2:length(data$dates)] - 1)

  api <- outbreaker2::get_cpp_api()

  ## Define default settings for o2mod.transphylo
  def_config <- outbreaker2::create_config(init_tree = init_tree,
                                           init_t_inf = init_t_inf,
                                           init_kappa = 1,
                                           init_pi = 1,
                                           find_import = FALSE,
                                           n_iter = 1000,
                                           sample_every = 1,
                                           move_mu = FALSE,
                                           move_pi = FALSE,
                                           move_eps = FALSE,
                                           move_lambda = FALSE,
                                           move_alpha = FALSE,
                                           move_swap_cases = TRUE,
                                           move_t_inf = TRUE,
                                           move_kappa = FALSE)
  def_priors <- outbreaker2::custom_priors()
  def_likelihoods <- outbreaker2::custom_likelihoods(genetic = lik_TransPhylo)
  def_moves <- outbreaker2::custom_moves(t_inf = new_move_tinf)

  ## Function to replace default settings with custom settings, if specified
  replace_def <- function(default, custom) {
    for(i in names(custom)) {
      default[[i]] <- custom[[i]]
    }
    return(default)
  }

  if(!is.null(c(config, priors, likelihoods, moves))) {
    warning("Modifying default settings may lead to unexpected behaviour; use with caution")
  }

  config <- replace_def(def_config, config)
  priors <- replace_def(def_priors, priors)
  likelihoods <- replace_def(def_likelihoods, likelihoods)
  moves <- replace_def(def_moves, moves)

  res <- outbreaker2::outbreaker(data = data ,
                                 config = config,
                                 priors = priors,
                                 likelihoods = likelihoods,
                                 moves = moves)

  return(res)

}
