library(tidyverse)
library(magrittr)
library(mulea)
source("source/functions.R")

#xxxxxxxxxxxxxxxxxxxxxxx
# AIM: calculating mean and SD of overlapping genes in all GMT files ---------------
#xxxxxxxxxxxxxxxxxxxxxxx
# to chose a highly overlapping gene set for eFDR and BH testing

#xxxxxxxxxxxxxxxxxxxxxxxx
# Read all GMT files ------------------------------------------------------
#xxxxxxxxxxxxxxxxxxxxxxxx
# Set the location of the folders containing the GMT files of each taxon
GMT_dir = "/home/barizona/Eszter/Kutatas/MulEA/GMT_files/GMT_files_for_mulea-main/GMT_files"
GMT_dir = "../../GMT_files/GMT_files_for_mulea-main/GMT_files"

# names of all GMT files
GMT_file_names <- list.files(path = GMT_dir, recursive = TRUE, full.names = TRUE)

# delete file names starting with "Genomic_location_Ensembl"
GMT_file_names <- GMT_file_names[!grepl("Genomic_location_Ensembl", GMT_file_names)]

# use only the uniprot ID files
GMT_file_names <- GMT_file_names[grepl("niprot", GMT_file_names)]

length(GMT_file_names)
# [1] 138

#xxxxxxxxxxxxxxxxxxxxxxxx
# Calculating the stats for all GMTs --------------------------------------
#xxxxxxxxxxxxxxxxxxxxxxxx

# create an empty matrix to store the results
GMT_statistics_matrix <- matrix(nrow = length(GMT_file_names), ncol = 7)
colnames(GMT_statistics_matrix) <- c("nr_of_terms",
                                     "mean_nr_of_elements",
                                     "sd_nr_of_elements",
                                     "median_nr_of_elements",
                                     "mean_nr_of_common_elements_between_2_terms",
                                     "sd_nr_of_common_elements_between_2_terms",
                                     "median_nr_of_common_elements_between_2_terms")
# loop through all GMT files
for (i in 1:length(GMT_file_names)) {
  # read the GMT file
  GMT <- read_gmt(file = GMT_file_names[i]) %>%
    filter_ontology(., min_nr_of_elements = 3, max_nr_of_elements = 400)

  # run the function only if the GMT remained big enough
  if (nrow(GMT) < 25) {
    print(paste("skip ", i))
    next
  }

  # calculate the statistics & store the results in the matrix
  GMT_statistics_matrix[i, ] <- GMT_statistics_func(GMT_variable = GMT)

  # print every i
  print(i)
}

rm(i, GMT)

# add row names to GMT_statistics_matrix
GMT_statistics_matrix %<>%
  data.frame() %>%
  tibble() %>%
  # names of all GMT files, cut the GMT_dir from the file names
  mutate(GMT_file_names = str_remove(GMT_file_names, paste0(GMT_dir, "/"))) %>%
  # relocate the GMT_file_names to the beginning
  relocate(GMT_file_names, .before = nr_of_terms)

# top 1 GMT file with the highest mean number of elements
GMT_statistics_matrix %>%
  arrange(desc(mean_nr_of_elements)) %>%
  head(1) %>%
  select(GMT_file_names, mean_nr_of_common_elements_between_2_terms)
# GMT_file_names                                                                                  mean_nr_of_common_el…¹
# <chr>                                                                                                            <dbl>
# 1 Saccharomyces_cerevisiae_4932/Transcription_factor_Yeastract_Saccharomyces_cerevisiae_UniprotI…                   16.2

# mean of mean of all ontologies
mean(GMT_statistics_matrix$mean_nr_of_common_elements_between_2_terms,
       na.rm = TRUE)
# 1.026056
