#' infer the time calibrated phylogeny associated with the provided alignment.
#' This function uses the R package babette to infer the phylogeny using BEAST2.
#' @param alignment Phydat object containing the focal alignment
#' @param treatment_name string to be appended to BEAST files
#' @param tree_prior tree prior used, default = birth-death prior
#' @param clock_prior clock prior used, default = strict clock
#' @param mcmc_seed seed of the mcmc chain, default is the system time
#' @param chain_length length of the mcmc chain, default is 1e7.
#' @param sample_interval interval of sampling, default is 5000
#' @param burnin burnin of posterior distribution
#' @param working_dir beast2 working dir
#' @param sub_rate substitution rate used to generate the original
#' alignment (if available), default is 1
#' @return list with all trees, and the consensus tree
#' @export
infer_phylogeny <- function(alignment,
                            treatment_name,
                            tree_prior = beautier::create_bd_tree_prior(),
                            clock_prior = beautier::create_strict_clock_model(),
                            mcmc_seed = NULL,
                            chain_length = 1e7,
                            sample_interval = 5000,
                            burnin = 0.1,
                            working_dir = NULL,
                            sub_rate = 1)  {

  if (is.null(working_dir)) working_dir <- getwd()

  if (!beastier::is_beast2_installed()) {
    stop("BEAST2 was not installed, please have a look at the package\n
         beastierinstall on GitHub: richelbilderbeek/beastierinstall")
  }

  temp_file_name <- "temp.fasta"
  phangorn::write.phyDat(alignment, file = temp_file_name, format = "fasta")

  if (is.null(mcmc_seed)) mcmc_seed <- round(as.numeric(Sys.time()))

  output_trees_filenames <- paste0(treatment_name, ".trees")
  output_log_filename <- paste0(treatment_name, ".log")

  inf_model <- beautier::create_inference_model(
    site_model = beautier::create_jc69_site_model(),
    clock_model = clock_prior,
    tree_prior = tree_prior,
    mcmc = beautier::create_mcmc(chain_length = chain_length,
                                 treelog =
                                   beautier::create_treelog(
                                              filename = output_trees_filenames,
                                              log_every = sample_interval),
                                 tracelog = beautier::create_tracelog(
                                              filename = output_log_filename,
                                              log_every = sample_interval)
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

  esses <- tracerer::calc_esses(beast_log, sample_interval = sample_interval)
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
