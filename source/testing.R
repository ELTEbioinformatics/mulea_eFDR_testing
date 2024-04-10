source("source/functions.R")

#xxxxxxxxxxxxxxxxxx
# AIM: running the tests with a GMT file ----------------------------------
#xxxxxxxxxxxxxxxxxx
# on the server:

# read & filter
GMT <- mulea::read_gmt(file = "input/Transcription_factor_Yeastract_Saccharomyces_cerevisiae_UniprotID.gmt") %>%
  filter_ontology(., min_nr_of_elements = 3, max_nr_of_elements = 400)

#xxxxxxxxxxxxxxxxxxxxxxxxxxxx
# Tests calculations ---------------------------------------------------------
#xxxxxxxxxxxxxxxxxxxxxxxxxxxx
sim_mult_tests_res <- simulateMultipleTestsWithRatioParam(
  input_gmt_filtered = GMT,
  noise_ratio_range = seq(0.0, 0.8, 0.1),
  over_repr_ratio = 0.85,
  number_of_tests = 1000,
  nthreads = 16)

saveRDS(sim_mult_tests_res, "output/sim_mult_tests_res_more_noise.Rds")

sim_mult_tests_res_sum <- getMultipleTestsSummaryAcrossCutOff(
  tests_res = sim_mult_tests_res,
  cut_off_range = seq(0, 1, 0.1))

saveRDS(sim_mult_tests_res_sum, "output/sim_mult_tests_res_sum_more_noise.Rds")

sim_mult_tests_res_roc <- getSummaryToRoc(sim_mult_tests_res)

saveRDS(sim_mult_tests_res_to_roc, "output/sim_mult_tests_res_to_roc_more_noise.Rds")

save.image()

q()





