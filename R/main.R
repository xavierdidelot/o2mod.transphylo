#' Run outbreaker2 using the TransPhylo model
#' @param data Data
#' @param n_iter Number of iterations
#' @param sample_every Thining interval
#' @return results of the run as an outbreaker2 output object
#' @export
o2mod.transphylo <- function(data,
                             n_iter = 1000,
                             sample_every = 1) {
  lik_TransPhylo <- function(data, param) {
    ttree <- list(
      ttree = cbind(param$t_inf, data$dates, param$alpha),
      nam = data$ptree$nam
    )
    ttree$ttree[which(is.na(ttree$ttree[, 3])), 3] <- 0
    txt <- utils::capture.output(ctree <- combine(ttree, data$ptree))
    if (length(txt) > 0)
      return(-Inf)
    return(TransPhylo::probPTreeGivenTTree(ctree, neg = 365 * 0.25))
  }

  deduceAlpha <- function(ptree, dates, alpha, t_inf) {
    ttree <- list(ttree = cbind(t_inf, dates, alpha), nam = ptree$nam)
    ttree$ttree[which(is.na(ttree$ttree[, 3])), 3] <- 0
    txt <- utils::capture.output(ctree <- combine(ttree, ptree))
    if (length(txt) > 0 && nchar(txt) > 9)
      return(alpha)
    alpha2 = TransPhylo::extractTTree(ctree)$ttree[, 3]
    return(replace(alpha2, alpha2 == 0, NA))
  }

  api <- outbreaker2::get_cpp_api()
  new_model <-
    outbreaker2::custom_likelihoods(genetic = lik_TransPhylo)
  new_move_tinf <-
    function(param, data, list_custom_ll = new_model) {
      for (i in 1:data$N) {
        current_ll <- api$cpp_ll_all(data, param, i = NULL, list_custom_ll)
        modif <- sample(c(-100:-1, 1:100), 1)
        param$t_inf[i] <- param$t_inf[i] + modif
        current_alpha <- param$alpha
        param$alpha <-
          deduceAlpha(data$ptree, data$dates, param$alpha, param$t_inf)
        new_ll <-
          api$cpp_ll_all(data, param, i = NULL, list_custom_ll)
        if (log(stats::runif(1)) > (new_ll - current_ll)) {
          param$t_inf[i] <- param$t_inf[i] - modif
          param$alpha <- current_alpha
        }
      }
      return(param)
    }
  new_moves <- outbreaker2::custom_moves(t_inf = new_move_tinf)

  init_tree <- c(NA, rep(1, length(data$dates) - 1))
  init_t_inf <-
    c(min(data$ptree$ptree[, 1]) - 1, data$dates[2:length(data$dates)] - 1)

  conf <-
    outbreaker2::create_config(
      init_tree = init_tree,
      init_t_inf = init_t_inf,
      init_kappa = 1,
      init_pi = 1,
      find_import = FALSE,
      n_iter = n_iter,
      sample_every = sample_every,
      move_mu = FALSE,
      move_pi = FALSE,
      move_eps = FALSE,
      move_lambda = FALSE,
      move_alpha = FALSE,
      move_swap_cases = TRUE,
      move_t_inf = TRUE,
      move_kappa = FALSE
    )

  res <- outbreaker2::outbreaker(
    data = data,
    config = conf,
    likelihoods = new_model,
    moves = new_moves
  )
}
