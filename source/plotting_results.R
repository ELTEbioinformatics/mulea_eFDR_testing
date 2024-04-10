library(tidyverse)
library(magrittr)
# library(ggpointdensity) # geom_pointdensity
# library(ggnewscale) # new_scale_color
# library(ggpubr) # ggarrange
# library(patchwork) # plot_layout
library(ggblend) # blend

#xxxxxxxxxxxxxxxxxxxxxx
# Tidy the tables ----------------------------------------------------
#xxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxx
## Summary table ----
#xxxxxxxxx
sim_mult_tests_res_sum %<>%
  # change names in method
  mutate(method = case_match(method,
                             "p_value" ~ "Uncorrected p-value",
                             "adjusted_p_value" ~ "Benjamini-Hochbeg FDR",
                             .default = method)) %>%
  # factorize method
  mutate(method = factor(method,
                         levels = c("Uncorrected p-value",
                                    "Benjamini-Hochbeg FDR",
                                    "eFDR")))

#xxxxxxxxx
## ROC table ----
#xxxxxxxxx
sim_mult_tests_res_roc %<>%
  # change names in method
  mutate(method = case_match(method,
                             "p_value" ~ "Uncorrected p-value",
                             "adjusted_p_value" ~ "Benjamini-Hochbeg FDR",
                             .default = method)) %>%
  # factorize method
  mutate(method = factor(method,
                         levels = c("Uncorrected p-value",
                                    "Benjamini-Hochbeg FDR",
                                    "eFDR")))

#xxxxxxxxxxxxxxxxxxxxxxxxxxx
# TPR - FPR plot -----------------------------------------------------------
#xxxxxxxxxxxxxxxxxxxxxxxxxxx

sim_mult_tests_res_sum %>%
  filter(method != "Uncorrected p-value") %>%
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
# png: width 700 height 900
# pdf: width 7 height 9 inch, cairo

# the resolution of the png is low
# to convert the pdf to good resolution png:
# convert -density 300 -trim TFLink_final_figure_v2.pdf -quality 100 TFLink_final_figure_v2b.png


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
rm(eFDR, BH, test, i)

t_test_pvals_TPRs


#xxxxxxxxxxxxxxxxxxxxx
# Precision - Recall plot -------------------------------------------------
#xxxxxxxxxxxxxxxxxxxxx

sim_mult_tests_res_sum %>%
  filter(method != "Uncorrected p-value") %>%
  # exclude noise 0.9
  filter(noise_ratio != 0.9) %>%
  # rewrite noise_ratio
  mutate(noise_ratio = paste0("Noise ratio: ", noise_ratio)) %>%
  ggplot(aes(x = Recall, y = Precision, color = method)) +
  geom_point(size = 0.1, alpha = 0.1)  * (blend("lighten") + blend("multiply", alpha = 0.1)) +
  # regression
  geom_line(stat = "smooth", method = "loess", formula = y ~ x, linewidth = 0.7) +
  scale_color_brewer(palette = "Set1", name = "p-value correction method") +
  facet_wrap( ~ noise_ratio) +
  theme_bw() +
  theme(legend.position = "bottom")


#xxxxxxxxxx
# ROC density curve plot -----------------------------------------------------
#xxxxxxxxxx
sim_mult_tests_res_roc %>%
  ggplot(aes(x = FPR,
             y = TPR,
             color = method)) +
  geom_point(alpha = 0.3) +
  stat_density_2d(geom = "point", aes(size = after_stat(density)), n = 20,
                  contour = FALSE, alpha = 0.3) +
  scale_radius(range = c(-0.5, 15)) +
  scale_color_brewer(palette = "Set1", name = "p-value correction method") +
  theme_bw() +
  theme(legend.position = "bottom")
