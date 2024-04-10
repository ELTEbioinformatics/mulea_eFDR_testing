source("source/functions.R")

#xxxxxxxxxxxxxxxxxx
# AIM: running the tests with a GMT file ----------------------------------
#xxxxxxxxxxxxxxxxxx
# on the server:


# read & filter
GMT <- mulea::filter_ontology(mulea::read_gmt(file = "input/test.GMT"),
                              min_nr_of_elements = 3, max_nr_of_elements = 400)

bg <- unique(unlist(GMT$list_of_values))

test_set <- sample(bg, 50, replace = FALSE)

#xxxxxxxxxxxxxxxxxxxxxxxxxxxx
# Tests calculations ---------------------------------------------------------
#xxxxxxxxxxxxxxxxxxxxxxxxxxxx
sim_mult_tests_res <- simulateMultipleTestsWithRatioParam(
  input_gmt_filtered = GMT,
  noise_ratio_range = seq(0.0, 0.8, 0.1),
  over_repr_ratio = 0.85,
  number_of_tests = 10,
  nthreads = 1)

sim_mult_tests_res_sum <- getMultipleTestsSummaryAcrossCutOff(
  tests_res = sim_mult_tests_res,
  cut_off_range = seq(0, 1, 0.1))

sim_mult_tests_res_roc <- getSummaryToRoc(sim_mult_tests_res)
