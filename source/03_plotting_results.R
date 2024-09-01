library(tidyverse)
library(magrittr)
library(ggblend) # blend

set.seed(1)

#xxxxxxxxxxxxxxxxxxxxxxxxxxx
# TPR - FPR plot: with a GMT file with highly overlapping entries --------------
#xxxxxxxxxxxxxxxxxxxxxxxxxxx

sim_mult_tests_res_sum_more_noise_GMT1 <- readRDS("output/sim_mult_tests_res_sum_more_noise_GMT1.Rds")

p <- sim_mult_tests_res_sum_more_noise_GMT1 %>%
  filter(method != "p_value") %>%
  # subsample percentage
  sample_frac(size = 0.5) %>%
  # change names in method
  mutate(method = case_match(method,
                             "adjusted_p_value" ~ "Benjamini-Hochbeg",
                             .default = method)) %>%
  # factorize method
  mutate(method = factor(method, levels = c("eFDR", "Benjamini-Hochbeg"))) %>%
  # rewrite noise_ratio
  mutate(noise_ratio = paste0("Noise ratio: ", noise_ratio)) %>%
  ggplot(aes(x = FPR, y = TPR, color = method)) +
  geom_point(size = 0.1, alpha = 0.1)  * (blend("lighten") + blend("multiply", alpha = 0.1)) +
  # regression
  geom_line(stat = "smooth", method = "loess", formula = y ~ x, linewidth = 0.7) +
  xlab("False positive rate") +
  ylab("True positive rate") +
  scale_color_brewer(palette = "Set1", name = "FDR method") +
  facet_wrap( ~ noise_ratio) +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.box="vertical")

ggsave("output/TPR_FPR_sim_mult_tests_res_sum_more_noise_GMT1.pdf", p,
       width = 4.6, height = 4.8, units = "in", device = cairo_pdf)

ggsave("output/TPR_FPR_sim_mult_tests_res_sum_more_noise_GMT1.png", p,
       width = 4.6, height = 4.8, units = "in", dpi = 300)

#xxxxxxxxxxxxxxxxxxxxxxxxxxx
# TPR - FPR plot: with a GMT file with lowly overlapping entries --------------
#xxxxxxxxxxxxxxxxxxxxxxxxxxx

sim_mult_tests_res_sum_more_noise_GMT2 <- readRDS("output/sim_mult_tests_res_sum_more_noise_GMT2.Rds")

p <- sim_mult_tests_res_sum_more_noise_GMT2 %>%
  filter(method != "p_value") %>%
  # subsample percentage
  sample_frac(size = 0.5) %>%
  # change names in method
  mutate(method = case_match(method,
                             "adjusted_p_value" ~ "Benjamini-Hochbeg",
                             .default = method)) %>%
  # factorize method
  mutate(method = factor(method, levels = c("eFDR", "Benjamini-Hochbeg"))) %>%
  # rewrite noise_ratio
  mutate(noise_ratio = paste0("Noise ratio: ", noise_ratio)) %>%
  ggplot(aes(x = FPR, y = TPR, color = method)) +
  geom_point(size = 0.1, alpha = 0.1)  * (blend("lighten") + blend("multiply", alpha = 0.1)) +
  # regression
  geom_line(stat = "smooth", method = "loess", formula = y ~ x, linewidth = 0.7) +
  xlab("False positive rate") +
  ylab("True positive rate") +
  scale_color_brewer(palette = "Set1", name = "FDR method") +
  facet_wrap( ~ noise_ratio) +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.box="vertical")

ggsave("output/TPR_FPR_sim_mult_tests_res_sum_more_noise_GMT2.pdf", p,
       width = 4.6, height = 4.8, units = "in", device = cairo_pdf)

ggsave("output/TPR_FPR_sim_mult_tests_res_sum_more_noise_GMT2.png", p,
       width = 4.6, height = 4.8, units = "in", dpi = 300)

