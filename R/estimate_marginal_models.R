#' infer the time calibrated phylogeny associated with the
#' @param fasta_filename file name of fasta file holding alignment for which the
#' marginal likelihood is to be estimated
#' @param rng_seed seed of pseudo-random number generator
#' @param verbose boolean indicating if verbose intermediate output is to be
#' generated
#' @return data frame
#' @export
estimate_marginal_models <- function(fasta_filename,
                                     rng_seed = 42,
                                     verbose = FALSE) {

  epsilon = 1e-12
  site_models = list(beautier::create_jc69_site_model())
  clock_models = beautier::create_clock_models()
  tree_priors = list(beautier::create_bd_tree_prior())

  if (rappdirs::app_dir()$os == "win") {
    stop("mcbette must run on Linux or Mac.\n", "\n", "It is not yet supported to call BEAST2 with packages installed\n",
         "in a scripted way")
  }
  if (!file.exists(fasta_filename)) {
    stop("'fasta_filename' must be the name of an existing FASTA file.\n",
         "File '", fasta_filename, "' not found")
  }
  beautier::check_site_models(site_models)
  beautier::check_clock_models(clock_models)
  beautier::check_tree_priors(tree_priors)
  if (!is.numeric(epsilon) || length(epsilon) != 1) {
    stop("'epsilon' must be one numerical value. Actual value(s): ",
         epsilon)
  }
  testit::assert(file.exists(fasta_filename))
  testit::assert(beastier::is_beast2_installed())
  testit::assert(mauricer::is_beast2_pkg_installed("NS"))
  n_rows <- length(site_models) * length(clock_models) * length(tree_priors)
  site_model_names <- rep(NA, n_rows)
  clock_model_names <- rep(NA, n_rows)
  tree_prior_names <- rep(NA, n_rows)
  marg_log_liks <- rep(NA, n_rows)
  marg_log_lik_sds <- rep(NA, n_rows)

  beast2_file_paths <- peregrine::create_pff_beast2_options()

  beast2_file_paths$output_log_filename <- "output_log.txt"
  beast2_file_paths$output_trees_filenames <- "output_trees.txt"
  beast2_file_paths$output_state_filename <- "output_state.txt"

  row_index <- 1
  for (site_model in site_models) {
    for (clock_model in clock_models) {
      for (tree_prior in tree_priors) {
        tryCatch({
          marg_lik <- babette::bbt_run(fasta_filename = fasta_filename,
                                       site_model = site_model,
                                       clock_model = clock_model,
                                       tree_prior = tree_prior,
                                       mcmc = beautier::create_mcmc_nested_sampling(epsilon = epsilon),
                                       beast2_path = beastier::get_default_beast2_bin_path(),
                                       beast2_output_log_filename = beast2_file_paths$output_log_filename,
                                       beast2_output_trees_filenames = beast2_file_paths$output_trees_filenames,
                                       beast2_output_state_filename = beast2_file_paths$output_state_filename,
                                       rng_seed = rng_seed,
                                       overwrite = TRUE)$ns
          marg_log_liks[row_index] <- marg_lik$marg_log_lik
          marg_log_lik_sds[row_index] <- marg_lik$marg_log_lik_sd
        }, error = function(msg) {
          if (verbose)
            print(msg)
        })
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

  calc_w <- function(v) {
    v <- v - max(v)
    a <- exp(-0.5 * v)
    b <- a / sum(a)
    return(b)

    #a <- exp(-0.5 *  Rmpfr::mpfr(marg_log_liks, 256))
    #b <- a / sum(a)
    #return(b)
  }

  weights <- as.numeric(calc_w(marg_log_liks) )

  df <- data.frame(site_model_name = site_model_names, clock_model_name = clock_model_names,
                   tree_prior_name = tree_prior_names, marg_log_lik = marg_log_liks,
                   marg_log_lik_sd = marg_log_lik_sds, weight = weights)
  df
}