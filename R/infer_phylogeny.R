#' infer the time calibrated phylogeny associated with the
#' @param alignment Phydat object containing the focal alignment
#' @param treatment_name string to be appended to BEAST files
#' @param tree_prior tree prior used
#' @param mcmc_seed seed of the mcmc chain, default is the system time
#' @param burnin burnin of posterior distribution
#' @param working_dir beast2 working dir
#' @param sub_rate substitution rate used to generate the original
#' alignment (if available), default is 1
#' @return list with all trees, and the consensus tree
#' @export
infer_phylogeny <- function(alignment,
                            treatment_name,
                            tree_prior = beautier::create_bd_tree_prior(),
                            mcmc_seed = NULL,
                            burnin = 0.1,
                            working_dir = NULL,
                            sub_rate = 1)  {

  if (is.null(working_dir)) working_dir <- getwd()

  temp_file_name <- "temp.fasta"
  phangorn::write.phyDat(alignment, file = temp_file_name, format = "fasta")

  if (is.null(mcmc_seed)) mcmc_seed <- round(as.numeric(Sys.time()))

  output_trees_filenames <- paste0(treatment_name, ".trees")
  output_log_filename <- paste0(treatment_name, ".log")

  inf_model <- beautier::create_inference_model(
    site_model = beautier::create_jc69_site_model(),
    clock_model = beautier::create_strict_clock_model(
      clock_rate_param = beautier::create_clock_rate_param(value = sub_rate),
      clock_rate_distr = beautier::create_gamma_distr()
    ),
    tree_prior = tree_prior,
    mcmc = beautier::create_mcmc(chain_length = 1e7,
                                 treelog = beautier::create_treelog(filename = output_trees_filenames,
                                                                    log_every = 5000),
                                 tracelog = beautier::create_tracelog(filename = output_log_filename,
                                                                      log_every = 5000)
    )
  )

  beast2_options_local <- beastier::create_beast2_options()

  posterior <- babette::bbt_run_from_model(
    fasta_filename = temp_file_name,
    inference_model = inf_model,
    beast2_options = beast2_options_local
  )

  file.remove(temp_file_name)

  beast_log <- tracerer::remove_burn_ins(
    posterior$estimates,
    burn_in_fraction = burnin
  )

  esses <- tracerer::calc_esses(beast_log, sample_interval = 5000)
  if (min(esses) < 200) {
    warning("WARNING! MCMC chain not converged!\n")
  }

  found_trees <- posterior$temp_trees

  remaining <- floor(burnin * length(found_trees))
  cat("remaining", remaining, "\n")
  cat("total_trees", length(found_trees), "\n")
  found_trees <- found_trees[remaining:length(found_trees)]

  consensus_tree <-  phangorn::maxCladeCred(found_trees)

  file.remove(output_log_filename)
  file.remove(output_trees_filenames)
  if (file.exists(beast2_options_local$input_filename)) {
    file.remove(beast2_options_local$input_filename)
  }
  if (file.exists(beast2_options_local$output_state_filename)) {
    file.remove(beast2_options_local$output_state_filename)
  }

  output <- list("all_trees" = found_trees,
                 "mcc_tree"  = consensus_tree)

  return(output)
}
