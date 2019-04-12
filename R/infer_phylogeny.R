#' infer the time calibrated phylogeny associated with the
#' @param alignment Phydat object containing the focal alignment
#' @param treatment_name string to be appended to BEAST files
#' @param burnin burnin of posterior distribution
#' @param chain_length length of the MCMC chain
#' @return list with all trees, and the consensus tree
#' @export
infer_phylogeny <- function(alignment,
                            treatment_name,
                            burnin,
                            chain_length)  {

  temp_file_name <- "temp.fasta"
  phangorn::write.phyDat(alignment, file = temp_file_name, format = "fasta")

  posterior <- babette::bbt_run_from_model(
    temp_file_name,
    inference_model = beautier::create_inference_model(
      mcmc = beautier::create_mcmc(chain_length = chain_length,
                                   store_every = 5000)
    ),
    beast2_options = beastier::create_beast2_options(
      overwrite = TRUE,
      output_trees_filenames = paste0(treatment_name, ".trees"),
      output_log_filename = paste0(treatment_name, ".log")
    )
  )

  file.remove(temp_file_name)

  beast_log <- tracerer::remove_burn_ins(
    posterior$estimates,
    burn_in_fraction = burnin
  )

  esses <- tracerer::calc_esses(beast_log, sample_interval = 5000)
  if (min(esses) < 200) {
    cat("WARNING! MCMC chain not converged!\n")
  }

  found_trees <- tracerer::parse_beast_trees(paste0(treatment_name, ".trees"))
  remaining <- floor(burnin * length(found_trees))
  found_trees <- found_trees[remaining:length(found_trees)]

  consensus_tree <-  phangorn::maxCladeCred(found_trees)

  output <- list("all_trees" = found_trees,
                 "mcc_tree" = consensus_tree)

  return(output)
}
