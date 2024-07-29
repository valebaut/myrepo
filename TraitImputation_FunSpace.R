library(funspace)
library(tidyverse)
library(dplyr)
library(readxl)
library(phytools)
library(ape)
library(phangorn)

# Read and process data
data.aux <- read_excel('C:/Users/valeb/OneDrive - Florida International University/data_species/CoralFunctionalTraits_PR.xlsx') %>%
  filter(trait_name %in% c("Skeletal density", "Growth rate", "Corallite diameter", "Colony maximum diameter", "Symbiodinium density")) %>%
  mutate(
    specie_name = str_trim(specie_name),
    trait_name = str_trim(trait_name) %>% tolower(),
    value = as.numeric(value)
  ) %>%
  mutate(
    standard_unit = case_when(
      trait_name == "skeletal density" & standard_unit %in% c("gCaCO3 cm^-3", "g CaCO3 cm^-3") ~ "g cm^-3",
      trait_name == "growth rate" & standard_unit == "mm d^-1" ~ "mm yr^-1",
      trait_name == "growth rate" & standard_unit == "mm month^-1" ~ "mm yr^-1",
      trait_name == "growth rate" & standard_unit == "mm yr" ~ "mm yr^-1",
      trait_name == "growth rate" & standard_unit == "cm yr" ~ "mm yr^-1",
      trait_name == "colony maximum diameter" & standard_unit == "mm" ~ "cm",
      trait_name == "symbiodinium density" & standard_unit == "cells x 10^6 cm^-2" ~ "cm^-2",
      TRUE ~ standard_unit
    ),
    value = case_when(
      trait_name == "growth rate" & standard_unit %in% c("mm d^-1", "mm month^-1") ~ value * ifelse(standard_unit == "mm d^-1", 365, 12),
      trait_name == "growth rate" & standard_unit == "cm yr" ~ value * 10,
      trait_name == "colony maximum diameter" & standard_unit == "mm" ~ value / 10,
      trait_name == "symbiodinium density" & standard_unit == "cells x 10^6 cm^-2" ~ value / 1e6,
      TRUE ~ value
    )
  ) %>%
  filter(!(standard_unit %in% c("cm^2/month", "g CaCO3 cm^-2", "% wt. gain month-1", "Âµm cm^-2 h^-1"))) %>%
  select(specie_name, trait_name, value) %>%
  group_by(specie_name, trait_name) %>%
  summarize(mean = mean(value, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = trait_name, values_from = mean) %>%
  mutate(specie_name = gsub(" ", "_", specie_name)) |>  
  mutate(across(2:6, ~ scale(log(.)), .names = "scaled_{col}")) |>  ##Scale and log transform 
  select(specie_name, starts_with("scaled_"))
  

##Adding phylogeny tree from Huang and Roy (2015) adapted by Carturan et al., (2021). File in Carturan Github.  

multi.phylo.tree<-read.tree("C:/Users/valeb/OneDrive - Florida International University/GitHubRep/PhylogeneticTree/supertree_names_validated_1000_798sp.tre")

print(class(multi.phylo.tree))  # Should be "multiPhylo"
print(length(multi.phylo.tree))  # Number of trees in the list


# Define the species list from your data
data_species <- unique(data.aux$specie_name)

# Function to find the best-matching tree
find_best_matching_tree <- function(trees, species_list) {
  best_tree <- NULL
  max_matches <- 0
  
  for (i in seq_along(trees)) {
    tree_species <- trees[[i]]$tip.label
    matches <- sum(species_list %in% tree_species)
    
    if (matches > max_matches) {
      best_tree <- trees[[i]]
      max_matches <- matches
    }
  }
  
  return(best_tree)
}

# Find the best-matching tree
phylo.tree <- find_best_matching_tree(multi.phylo.tree, data_species)
##check structure 
print(length(phylo.tree))
class(phylo.tree)

##structuring the trait data this way to remove species as a "column" 
data.aux <- column_to_rownames(data.aux, var = "specie_name")

##funspace package function combining phylogenetic tree and trait data with missing values 
traitinput<- impute(traits= data.aux,phylo= phylo.tree, addingSpecies=TRUE)

# Convert each element to a data frame
df_imputted <- as.data.frame(traitinput[[1]])
df_Original_Eigenvectors <- as.data.frame(traitinput[[2]])





