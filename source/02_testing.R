library(tidyverse)
library(magrittr)
library(mulea)
source("source/functions.R")

#xxxxxxxxxxxxxxxxxx
# AIM: running the tests with a GMT file with highly overlapping entries ------
#xxxxxxxxxxxxxxxxxx

# read & filter
GMT <- mulea::filter_ontology(
  mulea::read_gmt(file = "input/Transcription_factor_Yeastract_Saccharomyces_cerevisiae_UniprotID.gmt"),
  min_nr_of_elements = 3, max_nr_of_elements = 400)

# Tests calculations
sim_mult_tests_res_more_noise_GMT1 <- simulateMultipleTestsWithRatioParam(
  input_gmt_filtered = GMT,
  noise_ratio_range = c(0, 0.1, 0.2, 0.3),
  over_repr_ratio = 0.85,
  number_of_tests = 1000,
  nthreads = 4)

saveRDS(sim_mult_tests_res_more_noise_GMT1,
        "output/sim_mult_tests_res_more_noise_GMT1.Rds")

sim_mult_tests_res_sum_more_noise_GMT1 <- getMultipleTestsSummaryAcrossCutOff(
  tests_res = sim_mult_tests_res_more_noise_GMT1,
  cut_off_range = seq(0, 1, 0.1))

saveRDS(sim_mult_tests_res_sum_more_noise_GMT1,
        "output/sim_mult_tests_res_sum_more_noise_GMT1.Rds")

#xxxxxxxxxxxxxxxxxx
# AIM: running the tests with a GMT file with lowly overlapping entries ------
#xxxxxxxxxxxxxxxxxx

# read & filter
GMT <- mulea::filter_ontology(
  mulea::read_gmt(file = "input/GMT_files/Saccharomyces_cerevisiae_4932/Transcription_factor_TFLink_Saccharomyces_cerevisiae_SS_UniprotID.gmt"),
  min_nr_of_elements = 3, max_nr_of_elements = 400)

# Tests calculations
sim_mult_tests_res_more_noise_GMT2 <- simulateMultipleTestsWithRatioParam(
  input_gmt_filtered = GMT,
  noise_ratio_range = c(0, 0.1, 0.2, 0.3),
  over_repr_ratio = 0.85,
  number_of_tests = 1000,
  nthreads = 4)

saveRDS(sim_mult_tests_res_more_noise_GMT2,
        "output/sim_mult_tests_res_more_noise_GMT2.Rds")

sim_mult_tests_res_sum_more_noise_GMT2 <- getMultipleTestsSummaryAcrossCutOff(
  tests_res = sim_mult_tests_res_more_noise_GMT2,
  cut_off_range = seq(0, 1, 0.1))

saveRDS(sim_mult_tests_res_sum_more_noise_GMT2,
        "output/sim_mult_tests_res_sum_more_noise_GMT2.Rds")
