#xxxxxxxxxxxxxxxxxxxxxxxx
# calculating mean and SD of overlapping genes in GMT files ---------------
#xxxxxxxxxxxxxxxxxxxxxxxx

GMT_statistics_func <- function(GMT_variable) {
  # list of vectors containing elements
  GMT_element_list <- GMT_variable %>%
    select(list_of_values) %>%
    pull()

  # Lengths
  lengths <- sapply(GMT_element_list, length)

  # Intersections
  # all 2 length combinations of vectors
  combinations <- combn(length(GMT_element_list), 2)

  # Calculate the overlapping elements
  pairwise_intersections <- apply(combinations, 2, function(indices) {
    intersect(GMT_element_list[[indices[1]]], GMT_element_list[[indices[2]]])
  })

  # To get the lengths of the intersections
  pairwise_intersection_lengths <- sapply(pairwise_intersections, length)

  # Output
  # create a vector containing the statistics
  results <- c(
    nr_of_terms = length(GMT_element_list),
    mean_nr_of_elements = mean(lengths),
    sd_nr_of_elements = sd(lengths),
    median_nr_of_elements = median(lengths),
    mean_nr_of_common_elements_between_2_terms = mean(pairwise_intersection_lengths),
    sd_nr_of_common_elements_between_2_terms = sd(pairwise_intersection_lengths),
    median_nr_of_common_elements_between_2_terms = median(pairwise_intersection_lengths)
  )
  return(results)
}
# GMT <- read_gmt(file = "input/TFLink_Homo_sapiens_interactions_SS_GMT_proteinName_v1.0.gmt") %>%
#   filter_ontology(., min_nr_of_elements = 5, max_nr_of_elements = 400)
#
# GMT_statistics_func(GMT_variable = GMT)
# nr_of_terms
# 331.000000
# mean_nr_of_elements
# 34.190332
# sd_nr_of_elements
# 52.086089
# median_nr_of_elements
# 15.000000
# mean_nr_of_common_elements_between_2_terms
# 1.024389
# sd_nr_of_common_elements_between_2_terms
# 2.777999
# median_nr_of_common_elements_between_2_terms
# 0.000000

#xxxxxxxxxxxxxxxxxxxxx
# getMultipleTestsSummaryAcrossCutOff -------------------------------------
#xxxxxxxxxxxxxxxxxxxxx

getMultipleTestsSummaryAcrossCutOff <- function(tests_res,
                                                cut_off_range = seq(0, 1, 0.1)) {
  tests_res_sum <- NULL
  for (cut_off in cut_off_range) {
    tests_res_sum_p <- getMultipleTestsSummary(tests_res = tests_res,
                                               comparison_col_name = 'p_value',
                                               labels = list('method' = 'p', 'cut_off' = cut_off),
                                               cut_off = cut_off)
    tests_res_sum_bh <- getMultipleTestsSummary(tests_res = tests_res,
                                                comparison_col_name = 'adjustedPValue',
                                                labels = list('method' = 'bh', 'cut_off' = cut_off),
                                                cut_off = cut_off)
    tests_res_sum_pt <- getMultipleTestsSummary(tests_res = tests_res,
                                                comparison_col_name = 'eFDR',
                                                labels = list('method' = 'pt', 'cut_off' = cut_off),
                                                cut_off = cut_off)
    tests_res_sum <- rbind(tests_res_sum, tests_res_sum_p,
                           tests_res_sum_pt, tests_res_sum_bh)}
  return(tests_res_sum)
}


library(dplyr)
# PUBLIC API
#' @description
#' \code{decorateGmtByUnderOvenAndNoise}
#'
#' \code{decorateGmtByUnderOvenAndNoise} decorates GO with labels
#' (over, under, noise) per term.
#'
#' @param input_gmt input dataframe, read from gmt file.
#' @param number_of_over_representation_groups set the number of groups
#' which will be chosen to over representation.
#' @param number_of_under_representation_groups set the number of groups
#' which will be chosen to under representation.
#' @return Return data frame with model from specific location.
#' @noRd
decorateGmtByUnderOvenAndNoise <- function(input_gmt,
                                           number_of_over_representation_groups = 1,
                                           number_of_under_representation_groups = 0) {
  # Initialize all by noise labels.
  sample_label <- rep("noise", length(input_gmt$ontology_id))
  gmt_for_generator <- data.frame(input_gmt,
                                  "sample_label" = sample_label)
  # Choose and label terms for over and under representation.
  go_size <- length(gmt_for_generator$list_of_values)
  size_of_over_under_repr <- number_of_over_representation_groups + number_of_under_representation_groups
  go_change_repr <- sample(seq_len(go_size),
                           size_of_over_under_repr,
                           replace = FALSE)
  over_under_label <- c(rep("over", number_of_over_representation_groups),
                        rep("under", number_of_under_representation_groups))
  terms_to_manipulation <- data.frame('term_id' = go_change_repr,
                                      "over_under_label" = over_under_label)
  for (i in seq_along(terms_to_manipulation$term_id)) {
    term_row <- terms_to_manipulation[i, ]
    gmt_for_generator[term_row$term_id, ]$sample_label <- term_row$over_under_label}
  return(gmt_for_generator)
}


# PUBLIC API
#' @description
#' \code{generateInputSamples} Generates artificial GO with specific terms
#' (under or over represented).
#'
#' @param input_gmt input dataframe, read from gmt file.
#' @param noise_ratio level of noise in the sample, value from 0 to 1.
#' @param group_under_over_representation_ratio ratio of over represented group.
#' @param number_of_over_representation_groups number of groups chosen to over
#' representation.
#' @param number_of_under_representation_groups number of groups chosen to under
#' representation.
#' @return Return data frame with model from specific location.
#' @noRd
generateInputSamples <- function(input_gmt_decorated,
                                 noise_ratio = 0.2,
                                 over_repr_ratio = 0.85,
                                 under_repr_ratio = 0.05,
                                 rand_from_unique = TRUE,
                                 number_of_samples = 1) {
  all_genes_in_ontology <- NULL
  all_genes_in_enrichment <- NULL
  if (rand_from_unique) { all_genes_in_ontology <- unique(unlist(
    input_gmt_decorated$list_of_values))
  all_genes_in_enrichment <- unique(unlist(
    input_gmt_decorated[
      input_gmt_decorated$sample_label == "over", ]$list_of_values))
  } else {
    all_genes_in_ontology <- unlist(input_gmt_decorated$list_of_values)
    all_genes_in_enrichment <- unlist(input_gmt_decorated[
      input_gmt_decorated$sample_label == "over", ]$list_of_values)}
  size_of_ontology <- length(all_genes_in_ontology)
  size_of_noise <- ceiling(size_of_ontology * noise_ratio)
  size_of_enrichment <- ceiling(length(
    all_genes_in_enrichment) * over_repr_ratio)
  samples <- vector("list", number_of_samples)
  for (i in seq_along(samples)) {
    sample_noise <- all_genes_in_ontology[sample(seq_along(
      all_genes_in_ontology), size_of_noise, replace = FALSE)]
    sample_enrichment <- all_genes_in_enrichment[
      sample(seq_along(all_genes_in_enrichment), size_of_enrichment,
             replace = FALSE)]
    samples[[i]] <- unique(c(sample_noise, sample_enrichment))
    }
  return(samples)
}

# PUBLIC API
#' @description
#' \code{getSummaryToRoc}
#'
#' \code{getSummaryToRoc} generate artificial GO with specific terms under or
#' over represented.
#'
#' @param tests_res list of multiple tests results.
#' @return Return data frame which is the base to count ROC.
#' @noRd
#' @importFrom magrittr %>%
#' @importFrom plyr .
#' @importFrom rlang .data
getSummaryToRoc <- function(tests_res,
                            cut_off_resolution = 0.01,
                            methods_names = c("p_value", "adjusted_p_value", "eFDR")) {
  number_of_tests <- length(tests_res)
  data_to_roc <- data.frame("sample_label" = c(),
                            "p_value" = c(),
                            "adjusted_p_value" = c(),
                            "eFDR" = c())
  for (i in seq_len(number_of_tests)) {
    tests_res[[i]]$mulea_res[, c("p_value",
                                 "adjusted_p_value",
                                 "eFDR")]
    data_to_roc <- rbind(data_to_roc, data.frame("sample_label" = tests_res[[i]]$test_data[, c("sample_label")],
                                                 tests_res[[i]]$mulea_res[, c("p_value", "adjusted_p_value", "eFDR")]))
  }
  roc_stats <- tibble::tibble(TP_val = numeric(),
                              TN_val = numeric(),
                              FP_val = numeric(),
                              FN_val = numeric(),
                              TPR = numeric(),
                              FPR = numeric(),
                              sum_test = numeric(),
                              cut_off = numeric(),
                              method = character())
  for (method_name in methods_names) {
    for (cut_off in seq(0, 1, cut_off_resolution)) {
      sim_mult_tests_res_to_roc_summary <- data_to_roc %>%
        dplyr::mutate(PP = !!as.name(method_name) <= cut_off) %>%
        dplyr::mutate(TP = (.data$PP == TRUE & .data$sample_label == "over"),
                      TN = (.data$PP == FALSE & .data$sample_label != "over"),
                      FP = (.data$PP == TRUE & .data$sample_label != "over"),
                      FN = (.data$PP == FALSE & .data$sample_label == "over"))
      sim_sum <- sim_mult_tests_res_to_roc_summary %>%
        dplyr::summarise(TP_val = sum(.data$TP),
                         TN_val = sum(.data$TN),
                         FP_val = sum(.data$FP),
                         FN_val = sum(.data$FN))
      sim_sum_roc <- sim_sum %>%
        dplyr::mutate(TPR = .data$TP_val/(.data$TP_val + .data$FN_val),
                      FPR = .data$FP_val / (.data$FP_val + .data$TN_val),
                      sum_test = .data$TP_val + .data$TN_val + .data$FP_val + .data$FN_val,
                      cut_off = cut_off,
                      method = method_name)
      roc_stats <- roc_stats %>%
        tibble::add_row(sim_sum_roc)}}
  return(roc_stats)
}

# PUBLIC API
#' @description
#' \code{getMultipleTestsSummary}
#'
#' \code{getMultipleTestsSummary} generate artificial GO with specific terms
#' under or over represented.
#'
#' @param tests_res list of multiple tests results.
#' @param comparison_col_name column name which indicated data to compare on.
#' @param labels label datatable by additional columns with values.
#' @param cut_off threshold for value selected by comparison_col_name
#' @return Return data frame with FDR. TPRs per test.
#' @noRd
#' @importFrom magrittr %>%
#' @importFrom data.table :=
#' @import tictoc
getMultipleTestsSummary <- function(tests_res, comparison_col_name,
                                    labels = list(), cut_off = 0.05) {
  sumary_res <- prepareSummaryDf(tests_res)
  number_of_tests <- length(tests_res)
  for (i in seq_len(number_of_tests)) {
    total_population <- tests_res[[i]]$test_data$ontology_id
    total_population_size <- length(total_population)
    P <- tests_res[[i]]$test_data[
      tests_res[[i]]$test_data$sample_label == 'over', ]$ontology_id
    P_size <- length(P)
    N <- tests_res[[i]]$test_data[
      tests_res[[i]]$test_data$sample_label != 'over', ]$ontology_id
    N_size <- length(N)
    if (P_size + N_size != total_population_size) {
      warning("Not OK size of Actual in contingency table")}
    PP <- tests_res[[i]]$mulea_res[tests_res[[i]]$mulea_res[, comparison_col_name] <= cut_off,]$ontology_id
    PP_size <- length(PP)
    PN <- tests_res[[i]]$mulea_res[tests_res[[i]]$mulea_res[, comparison_col_name] > cut_off,]$ontology_id
    PN_size <- length(PN)
    if (PP_size + PN_size != total_population_size) {
      warning("Not OK size of Predicted in contingency table")}
    TP <- intersect(P, PP)
    TP_size <- length(TP)
    FP <- intersect(N, PP)
    FP_size <- length(FP)
    FN <- intersect(P, PN)
    FN_size <- length(FN)
    TN <- intersect(N, PN)
    TN_size <- length(TN)
    if (TP_size + FP_size + FN_size + TN_size != total_population_size) {
      warning("Not OK size of total contingency table")
      }
    over_repr_terms <- tests_res[[i]]$test_data[tests_res[[i]]$test_data$sample_label == "over", ]$ontology_id
    sumary_res_tmp <- getMetadata(i = i,
                                  TP = TP, TP_size = TP_size,
                                  FP = FP, FP_size = FP_size,
                                  FN = FN, FN_size = FN_size,
                                  TN = TN, TN_size = TN_size,
                                  over_repr_terms = over_repr_terms,
                                  tests_res = tests_res)
    sumary_res[i,] <- sumary_res_tmp
  }
  sumary_res <- getSummaryRes(sumary_res = sumary_res,
                              FP_size = FP_size,
                              TN_size = TN_size,
                              TP_size = TP_size,
                              FN_size = FN_size)
  sumary_res <- addLabelsToSummary(labels = labels,
                                   sumary_res = sumary_res)
  return(sumary_res)
}

prepareSummaryDf <- function(tests_res) {
  metadata_len <- length(tests_res[[1]]$metadata)
  sumary_res <- data.frame(matrix(ncol = 10 + metadata_len, nrow = 0))
  colnames(sumary_res) <- c("test_no",
                            "TP", "TP_size",
                            "FP", "FP_size",
                            "FN", "FN_size",
                            "TN", "TN_size",
                            "over_repr_terms",
                            names(tests_res[[1]]$metadata))
  return(sumary_res)
}

getSummaryRes <- function(sumary_res, FP_size, TN_size, TP_size, FN_size) {
  sumary_res <- tibble::tibble(sumary_res) %>%
    dplyr::mutate(FPR = FP_size / (FP_size + TN_size)) %>%
    dplyr::mutate(TPR = TP_size / (TP_size + FN_size)) %>%
    dplyr::mutate(FDR = FP_size / (TP_size + FP_size)) %>%
    dplyr::mutate(NPV = TN_size / (FN_size + TN_size)) %>%
    dplyr::mutate(Precision = TP_size / (TP_size + FP_size)) %>%
    dplyr::mutate(Recall = TP_size / (TP_size + FN_size))
  sumary_res
}

getMetadata <- function(i,
                        TP, TP_size,
                        FP, FP_size,
                        FN, FN_size,
                        TN, TN_size,
                        over_repr_terms, tests_res) {
  sumary_res_tmp <- data.frame("test_no" = i,
                               "TP" = I(list(TP)), "TP_size" = TP_size,
                               "FP" = I(list(FP)), "FP_size" = FP_size,
                               "FN" = I(list(FN)), "FN_size" = FN_size,
                               "TN" = I(list(TN)), "TN_size" = TN_size,
                               "over_repr_terms" = I(list(over_repr_terms)))
  for (metadata_entry in names(tests_res[[i]]$metadata)) {
    if ("input_select" == metadata_entry) {
      sumary_res_tmp <- cbind(sumary_res_tmp,
                              "input_select" = I(tests_res[[i]]$metadata[metadata_entry]))
    } else { sumary_res_tmp <- cbind(sumary_res_tmp,
                                     metadata_entry = as.character(tests_res[[i]]$metadata[metadata_entry]))
    }
  }
  sumary_res_tmp
}

addLabelsToSummary <- function(labels, sumary_res) {
  for (label_id in seq_along(labels)) {
    # IMPORTANT : Labels are as characters in datatable
    label_name <- as.character(names(labels)[[label_id]])
    label_value <- as.character(labels[[label_id]])
    sumary_res <- sumary_res %>%
      dplyr::mutate(!!label_name := label_value)}
  sumary_res
}

# PUBLIC API
#' @description
#' \code{getMultipleTestsSummaryAcrossCutOff}
#'
#' \code{getMultipleTestsSummaryAcrossCutOff} doing summary across cutoff range.
#'
#' @param tests_res list of multiple tests results.
#' @param cut_off_range threshold for value selected by comparison_col_name
#' @return Return data frame with FDR. TPRs per test.
#' @noRd
getMultipleTestsSummaryAcrossCutOff <- function(tests_res,
                                                cut_off_range = seq(0, 1, 0.1)) {
  tests_res_sum <- NULL
  for (cut_off in cut_off_range) {
    tests_res_sum_p <- getMultipleTestsSummary(tests_res = tests_res,
                                               comparison_col_name = 'p_value',
                                               labels = list('method' = 'p_value', 'cut_off' = cut_off),
                                               cut_off = cut_off)
    tests_res_sum_bh <- getMultipleTestsSummary(tests_res = tests_res,
                                                comparison_col_name = 'adjusted_p_value',
                                                labels = list('method' = 'adjusted_p_value', 'cut_off' = cut_off),
                                                cut_off = cut_off)
    tests_res_sum_pt <- getMultipleTestsSummary(tests_res = tests_res,
                                                comparison_col_name = 'eFDR',
                                                labels = list('method' = 'eFDR', 'cut_off' = cut_off),
                                                cut_off = cut_off)
    tests_res_sum <- rbind(tests_res_sum, tests_res_sum_p,
                           tests_res_sum_pt, tests_res_sum_bh)
    }
  return(tests_res_sum)
}


# PUBLIC API
#' @description
#' \code{simulateMultipleTests}
#'
#' \code{simulateMultipleTests} generate artificial GO with specific terms under
#' or over represented.
#'
#' @param input_gmt_filtered gmt data frame with ontology for tests.
#' @param number_of_tests number of tests to perform.
#' @param noise_ratio ratio of noise in data from [0,1] interval.
#' @param number_of_over_representation_groups number of terms to over
#' represent.
#' @param number_of_under_representation_groups number of terms to under
#' represent.
#' @return Return data frame with FDR. TPRs per test.
#' @noRd
simulateMultipleTests <- function(input_gmt_filtered,
                                  number_of_tests = 10,
                                  noise_ratio = 0.35,
                                  over_repr_ratio  = 0.85,
                                  number_of_over_representation_groups = ceiling(nrow(input_gmt_filtered) * 0.1),
                                  number_of_under_representation_groups = 0,
                                  number_of_steps = 5000, nthreads = 1) {
  tictoc::tic()
  number_of_samples <- 1
  tests_res <- vector("list", number_of_tests)
  for (i in seq_len(number_of_tests)) {
    input_gmt_decorated <- decorateGmtByUnderOvenAndNoise(input_gmt = input_gmt_filtered,
                                                          number_of_over_representation_groups = number_of_over_representation_groups,
                                                          number_of_under_representation_groups = number_of_under_representation_groups)
    samples <- generateInputSamples(input_gmt_decorated,
                                    noise_ratio = noise_ratio,
                                    over_repr_ratio = over_repr_ratio,
                                    number_of_samples = number_of_samples)
    if (length(samples) != 1) { warning("sample is not size 1") }
    input_select <- unlist(samples)
    # eFDR
    mulea_ora_model <- mulea::ora(gmt = input_gmt_filtered,
                                  element_names = input_select,
                                  p_value_adjustment_method = "eFDR",
                                  number_of_permutations = number_of_steps,
                                  nthreads = nthreads)
    mulea_ora_results <- mulea::run_test(mulea_ora_model)
    # BH
    mulea_ora_model_2 <- mulea::ora(gmt = input_gmt_filtered,
                                    element_names = input_select,
                                    p_value_adjustment_method = "BH",
                                    # number_of_permutations = number_of_steps,
                                    nthreads = 1)
    mulea_ora_results_2 <- mulea::run_test(mulea_ora_model_2)
    # Join of the 2 results
    mulea_ora_results <- merge(mulea_ora_results,
                             # Selecting columns from mulea_ora_results_2
                             mulea_ora_results_2[, c("ontology_id", "adjusted_p_value")],
                             by = "ontology_id", all.x = TRUE)

    tests_res[[i]]$mulea_res <- mulea_ora_results
    tests_res[[i]]$test_data <- input_gmt_decorated
    tests_res[[i]]$metadata <- list("noise_ratio" = noise_ratio,
                                    "number_of_tests" = number_of_tests,
                                    "over_repr_ratio" = over_repr_ratio,
                                    "number_of_over_representation_groups" = number_of_over_representation_groups,
                                    "input_select" = input_select)
    }
  tictoc::toc()
  return(tests_res)
}


# PUBLIC API
#' @description
#' \code{simulateMultipleTestsWithRatioParam}
#'
#' \code{simulateMultipleTestsWithRatioParam} generate artificial GO with
#' specific terms under or over represented.
#'
#' @param input_gmt_filtered gmt data frame with ontology for tests.
#' @param number_of_tests number of tests to perform.
#' @param noise_ratio_range range of ratios of noise in data from [0,1]
#' interval.
#' @param number_of_over_representation_groups number of terms to over
#' represent.
#' @param number_of_steps number of steps.
#' @return Return data frame with FDR. TPRs per test.
#' @noRd
simulateMultipleTestsWithRatioParam <- function(input_gmt_filtered,
                                                noise_ratio_range = c(0.2, 0.5),
                                                number_of_tests = 10,
                                                over_repr_ratio = 0.85,
                                                number_of_over_representation_groups = ceiling(nrow(input_gmt_filtered) * 0.2),
                                                number_of_steps = 5000,
                                                nthreads = 1) {
  tictoc::tic()
  sim_mult_tests <- list()
  for (noise_ratio in noise_ratio_range) {
    sim_mult_tests <- c(sim_mult_tests,
                        simulateMultipleTests(input_gmt_filtered = input_gmt_filtered,
                                              number_of_tests = number_of_tests,
                                              noise_ratio = noise_ratio,
                                              over_repr_ratio = over_repr_ratio,
                                              number_of_over_representation_groups = number_of_over_representation_groups,
                                              number_of_under_representation_groups = 0,
                                              number_of_steps = number_of_steps,
                                              nthreads = nthreads))
    }
  tictoc::toc()
  return(sim_mult_tests)
}
