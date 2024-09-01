library(tidyverse)

#xxxxxxxxxxxxxxxxxxxxx
# Testing the difference between the TPRs: highly overlapping entries ---------
#xxxxxxxxxxxxxxxxxxxxx

t_test_pvals_TPRs_more_noise_GMT1 <- c()

for(i in unique(sim_mult_tests_res_sum_more_noise_GMT1$noise_ratio)) {
  #xxxxx
  # True positive rates
  #xxxxx
  eFDR <- sim_mult_tests_res_sum_more_noise_GMT1 %>%
    filter(method == "eFDR") %>%
    filter(noise_ratio == i) %>%
    select(TPR) %>%
    unlist(use.names = FALSE)

  BH <- sim_mult_tests_res_sum_more_noise_GMT1 %>%
    filter(method == "adjusted_p_value") %>%
    filter(noise_ratio == i) %>%
    select(TPR) %>%
    unlist(use.names = FALSE)

  test <- t.test(eFDR, BH, alternative = "greater", paired = TRUE)

  t_test_pvals_TPRs_more_noise_GMT1 <- c(t_test_pvals_TPRs_more_noise_GMT1,
                                          test$p.value)
}

names(t_test_pvals_TPRs_more_noise_GMT1) <- unique(sim_mult_tests_res_sum_more_noise_GMT1$noise_ratio)

t_test_pvals_TPRs_more_noise_GMT1 %>%
  enframe(name = "noise", value = "t-test_p-val") %>%
  write_tsv("output/t_test_pvals_TPRs_more_noise_GMT1.tsv")

rm(eFDR, BH, test, i)

#xxxxxxxxxxxxxxxxxxxxx
# Testing the difference between the TPRs: lowly overlapping entries ---------
#xxxxxxxxxxxxxxxxxxxxx

t_test_pvals_TPRs_more_noise_GMT2 <- c()

for(i in unique(sim_mult_tests_res_sum_more_noise_GMT2$noise_ratio)) {
  #xxxxx
  # True positive rates
  #xxxxx
  eFDR <- sim_mult_tests_res_sum_more_noise_GMT2 %>%
    filter(method == "eFDR") %>%
    filter(noise_ratio == i) %>%
    select(TPR) %>%
    unlist(use.names = FALSE)

  BH <- sim_mult_tests_res_sum_more_noise_GMT2 %>%
    filter(method == "adjusted_p_value") %>%
    filter(noise_ratio == i) %>%
    select(TPR) %>%
    unlist(use.names = FALSE)

  test <- t.test(eFDR, BH, alternative = "greater", paired = TRUE)

  t_test_pvals_TPRs_more_noise_GMT2 <- c(t_test_pvals_TPRs_more_noise_GMT2,
                                         test$p.value)
}

names(t_test_pvals_TPRs_more_noise_GMT2) <- unique(sim_mult_tests_res_sum_more_noise_GMT2$noise_ratio)

t_test_pvals_TPRs_more_noise_GMT2 %>%
  enframe(name = "noise", value = "t-test_p-val") %>%
  write_tsv("output/t_test_pvals_TPRs_more_noise_GMT2.tsv")

rm(eFDR, BH, test, i)
