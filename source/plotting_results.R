library(tidyverse)
library(magrittr)
library(ggblend) # blend

set.seed(1)
#xxxxxxxxxxxxxxxxxxxxxx
# Tidy the tables ----------------------------------------------------
#xxxxxxxxxxxxxxxxxxxxxx

sim_mult_tests_res_sum <- readRDS("output/sim_mult_tests_res_sum_more_noise.Rds")

#xxxxxxxxxxxxxxxxxxxxxxxxxxx
# TPR - FPR plot 0 - 0.3 noise level -----------------------------------------------------------
#xxxxxxxxxxxxxxxxxxxxxxxxxxx

sim_mult_tests_res_sum %>%
  filter(method != "p_value" & noise_ratio <= 0.3) %>%
  # subsample percentage
  sample_frac(size = 0.5) %>%
  # change names in method
  mutate(method = case_match(method,
                             "adjusted_p_value" ~ "Benjamini-Hochbeg FDR",
                             .default = method)) %>%
  # factorize method
  mutate(method = factor(method, levels = c("eFDR", "Benjamini-Hochbeg FDR"))) %>%
  # rewrite noise_ratio
  mutate(noise_ratio = paste0("Noise ratio: ", noise_ratio)) %>%
  ggplot(aes(x = FPR, y = TPR, color = method)) +
  geom_point(size = 0.1, alpha = 0.1)  * (blend("lighten") + blend("multiply", alpha = 0.1)) +
  # regression
  geom_line(stat = "smooth", method = "loess", formula = y ~ x, linewidth = 0.7) +
  xlab("False positive rate") +
  ylab("True positive rate") +
  scale_color_brewer(palette = "Set1", name = "p-value correction method") +
  facet_wrap( ~ noise_ratio) +
  theme_bw() +
  theme(legend.position = "bottom")

# can be saved manually only
# pdf: width 4.6 height 4.8 inch, cairo

# the resolution of the png is low
# to convert the pdf to good resolution png:
# convert -density 300 -trim output/TPR_FPR_subsample_0.5_noise_0-0.3.pdf -quality 100 output/TPR_FPR_subsample_0.5_noise_0-0.3.png


#xxxxxxxxxxxxxxxxxxxxxxxxxxx
# TPR - FPR plot 0 - 0.9 noise level -----------------------------------------------------------
#xxxxxxxxxxxxxxxxxxxxxxxxxxx

sim_mult_tests_res_sum %>%
  filter(method != "p_value") %>%
  # subsample percentage
  sample_frac(size = 0.5) %>%
  # change names in method
  mutate(method = case_match(method,
                             "adjusted_p_value" ~ "Benjamini-Hochbeg FDR",
                             .default = method)) %>%
  # factorize method
  mutate(method = factor(method, levels = c("eFDR", "Benjamini-Hochbeg FDR"))) %>%
  # rewrite noise_ratio
  mutate(noise_ratio = paste0("Noise ratio: ", noise_ratio)) %>%
  ggplot(aes(x = FPR, y = TPR, color = method)) +
  geom_point(size = 0.1, alpha = 0.1)  * (blend("lighten") + blend("multiply", alpha = 0.1)) +
  # regression
  geom_line(stat = "smooth", method = "loess", formula = y ~ x, linewidth = 0.7) +
  xlab("False positive rate") +
  ylab("True positive rate") +
  scale_color_brewer(palette = "Set1", name = "p-value correction method") +
  facet_wrap( ~ noise_ratio) +
  theme_bw() +
  theme(legend.position = "bottom")

# can be saved manually only
# pdf: width 7 height 7.5 inch, cairo

# the resolution of the png is low
# to convert the pdf to good resolution png:
# convert -density 300 -trim TPR_FPR_subsample_0.5.pdf -quality 100 TPR_FPR_subsample_0.5.png


#xxxxxxxxxxxxxxxxxxxxx
## Testing if there is a significant difference between the TPRs --------
#xxxxxxxxxxxxxxxxxxxxx

t_test_pvals_TPRs <- c()

for(i in unique(sim_mult_tests_res_sum$noise_ratio)) {
  #xxxxx
  # True positive rates
  #xxxxx
  eFDR <- sim_mult_tests_res_sum %>%
    filter(method == "eFDR") %>%
    # exclude noise 0.9
    filter(noise_ratio == i) %>%
    select(TPR) %>%
    unlist(use.names = FALSE)

  BH <- sim_mult_tests_res_sum %>%
    filter(method == "Benjamini-Hochbeg FDR") %>%
    # exclude noise 0.9
    filter(noise_ratio == i) %>%
    select(TPR) %>%
    unlist(use.names = FALSE)

  test <- t.test(eFDR, BH, alternative = "greater", paired = TRUE)

  t_test_pvals_TPRs <- c(t_test_pvals_TPRs, test$p.value)
}

names(t_test_pvals_TPRs) <- unique(sim_mult_tests_res_sum$noise_ratio)
# p-value < 2.2e-16 for all
rm(eFDR, BH, test, i, sim_mult_tests_res_sum)

t_test_pvals_TPRs
#             0           0.1           0.2           0.3           0.4
# 9.017472e-213 1.997708e-219 2.137203e-222 3.791461e-226 5.159119e-231
#           0.5           0.6           0.7           0.8
# 6.564688e-232 1.055894e-236 3.918375e-246 1.091873e-295

