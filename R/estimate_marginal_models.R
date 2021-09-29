get_marg_lik <- function(fasta_filename,
                         site_model,
                         clock_model,
                         tree_prior,
                         rng_seed) {

  beast2_input_filename <- beastier::create_temp_input_filename()
  beast2_output_state_filename <- beastier::create_temp_state_filename()

  marg_lik <- babette::bbt_run(
    fasta_filename = fasta_filename,
    site_model = site_model,
    clock_model = clock_model,
    tree_prior = tree_prior,
    beast2_input_filename = beast2_input_filename,
    beast2_output_state_filename = beast2_output_state_filename,
    mcmc =
      beautier::create_ns_mcmc(chain_length = 1e9,
               store_every = 5000,
               tracelog = beautier::create_tracelog(filename = "marg.trace"),
               treelog = beautier::create_treelog(filename = "marg.trees")),
    beast2_path = beastier::get_default_beast2_bin_path(),
    rng_seed = rng_seed,
    overwrite = TRUE)$ns

  file.remove(beast2_output_state_filename)
  file.remove(beast2_input_filename)

  return(marg_lik)
}

#' estimate the marginal likelihood of the relaxed and strict clock model for
#' a provided alignment
#' @description estimate_marginal_models estimates the marginal likelihood of
#' both the strict and the relaxed clock model, given the JC69 substitution
#' model, using the NS package in BEAST, made available via the babette R
#' package. The NS package performs nested sampling, and uses an MCMC approach
#' to estimate the marginal likelihood. Sampling is performed until convergence
#' of the MCMC chain.
#' @param fasta_filename file name of fasta file holding alignment for which the
#' marginal likelihood is to be estimated
#' @param use_yule_prior by default, a birth-death prior is used as tree prior,
#' but if use_yule_prior is set to TRUE, a pure-birth prior will be used.
#' @param rng_seed seed of pseudo-random number generator
#' @param sub_rate substitution rate
#' @param verbose boolean indicating if verbose intermediate output is to be
#' generated
#' @return data frame with marginal likelihoods and relative weights per clock
#' model.
#' @export
estimate_marginal_models <- function(fasta_filename,
                                     use_yule_prior = FALSE,
                                     rng_seed = 42,
                                     sub_rate = 1.0,
                                     verbose = FALSE) {

  site_models <- list(beautier::create_jc69_site_model())
  clock_models <- list()
  clock_models[[1]] <- beautier::create_rln_clock_model(
                          mean_clock_rate = sub_rate,
                          mean_rate_prior_distr =
                                beautier::create_distr_log_normal())
  clock_models[[2]] <- beautier::create_strict_clock_model(
                          clock_rate_param =
                            beautier::create_clock_rate_param(value = sub_rate),
                          clock_rate_distr =
                            beautier::create_distr_log_normal())
  tree_priors <- list(beautier::create_bd_tree_prior())
  if (use_yule_prior) {
    tree_priors <- list(beautier::create_yule_tree_prior())
  }

  if (rappdirs::app_dir()$os == "win") {
    stop("mcbette must run on Linux or Mac.\n",
         "\n",
         "It is not yet supported to call BEAST2 with packages installed\n",
         "in a scripted way")
  }
  if (!file.exists(fasta_filename)) {
    stop("'fasta_filename' must be the name of an existing FASTA file.\n",
         "File '", fasta_filename, "' not found")
  }

  beautier::check_site_models(site_models)
  beautier::check_clock_models(clock_models)
  beautier::check_tree_priors(tree_priors)

  testit::assert(file.exists(fasta_filename))
  testit::assert(beastier::is_beast2_installed())
  testit::assert(mauricer::is_beast2_pkg_installed("NS"))

  n_rows <- length(site_models) *
            length(clock_models) *
            length(tree_priors)

  site_model_names  <- rep(NA, n_rows)
  clock_model_names <- rep(NA, n_rows)
  tree_prior_names  <- rep(NA, n_rows)
  marg_log_liks     <- rep(NA, n_rows)
  marg_log_lik_sds  <- rep(NA, n_rows)

  row_index <- 1

  # this code looks awful right now. TODO: clean up!

  for (site_model in site_models) {
    for (clock_model in clock_models) {
      for (tree_prior in tree_priors) {
        tryCatch({
          marg_lik <- get_marg_lik(fasta_filename,
                                   site_model,
                                   clock_model,
                                   tree_prior,
                                   rng_seed)
          marg_log_liks[row_index] <- marg_lik$marg_log_lik
          marg_log_lik_sds[row_index] <- marg_lik$marg_log_lik_sd},
          error = function(msg) {
            if (verbose)
              print(msg)
          }
        )
        site_model_names[row_index] <- site_model$name
        clock_model_names[row_index] <- clock_model$name
        tree_prior_names[row_index] <- tree_prior$name
        if (verbose == TRUE) {
          print(paste0("Log evidence for model ", row_index, "/", n_rows,
                       ": ", marg_log_liks[row_index]))
        }
        row_index <- row_index + 1
      }
    }
  }

  # following code was ripped from mcbette, but mcbette is not on CRAN yet.
  marg_liks <- exp(Rmpfr::mpfr(marg_log_liks, 256))
  weights <- rep(Rmpfr::mpfr(0, 256), length(marg_liks))
  if (sum(marg_liks) != Rmpfr::mpfr(0, 256)) {
    weights <- marg_liks / sum(marg_liks)
  }

  df <- data.frame(site_model_name = site_model_names,
                   clock_model_name = clock_model_names,
                   tree_prior_name = tree_prior_names,
                   marg_log_lik = marg_log_liks,
                   marg_log_lik_sd = marg_log_lik_sds,
                   weight = as.numeric(weights))

  if (file.exists("marg.trace")) {
    file.remove("marg.trace")
  }
  if (file.exists("marg.trees")) {
    file.remove("marg.trees")
  }
  if (file.exists("marg.posterior.trace")) {
    file.remove("marg.posterior.trace")
  }
  if (file.exists("marg.posterior.trees")) {
    file.remove("marg.posterior.trees")
  }

  return(df)
}
