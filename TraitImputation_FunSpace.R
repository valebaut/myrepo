library(funspace)
library(tidyverse)
library(dplyr)
library(readxl)
df_tr = read_excel('C:/Users/valeb/OneDrive - Florida International University/data_species/CoralFunctionalTraits_PR.xlsx') |>
  filter(trait_name %in% c("Skeletal density","Growth rate", "Corallite diameter", "Colony maximum diameter", "Symbiodinium density")) |> 
  mutate(specie_name = str_trim(specie_name),
         trait_name=str_trim(trait_name),
         trait_name=tolower(trait_name),
         value=as.numeric(value))  |> 
  mutate(standard_unit = ifelse(trait_name == "skeletal density" & standard_unit == "gCaCO3 cm^-3", "g cm^-3", standard_unit),
         standard_unit = ifelse(trait_name == "skeletal density" & standard_unit == "g CaCO3 cm^-3", "g cm^-3", standard_unit),
         # ^ changing unit only (value stay the same) 
         value = ifelse(trait_name == "growth rate" & standard_unit == "mm d^-1", value * 365, value),
         standard_unit = ifelse(trait_name == "growth rate" & standard_unit == "mm d^-1", "mm yr^-1", standard_unit),
         # ^ changing value and standard unit for mm/day  -> mm/year
         value = ifelse(trait_name == "growth rate" & standard_unit == "mm month^-1", value * 12, value),
         standard_unit = ifelse(trait_name == "growth rate" & standard_unit == "mm month^-1", "mm yr^-1", standard_unit),
         # ^ changing value and standard unit for mm/month -> mm/year
         standard_unit = ifelse(trait_name == "growth rate" & standard_unit == "mm yr", "mm yr^-1", standard_unit),
         # ^ changing unit only (value stay the same)
         value = ifelse(trait_name == "growth rate" & standard_unit == "cm yr", value * 10, value), 
         standard_unit = ifelse(trait_name == "growth rate" & standard_unit == "cm yr", "mm yr^-1", standard_unit),
         # ^ changing value and standard unit for cm/year -> mm/year
         value = ifelse(trait_name == "colony maximum diameter"  & standard_unit == "mm", value/10, value),
         standard_unit = ifelse(trait_name == "colony maximum diameter"  & standard_unit == "mm", "cm"  , standard_unit),
         # ^ changing value and standard unit for  "mm" -> "cm"
         # value = ifelse(trait_name == "calcification rate" & standard_unit == "mg cm^-2 d^-1", value*0.365, value),
         # standard_unit = ifelse(trait_name == "calcification rate" & standard_unit == "mg cm^-2 d^-1", "g cm^-2 yr^-1", standard_unit),
         # # ^ changing value and standard unit for mg/cm^2/day -> g/cm^2/yr 
         # value = ifelse(trait_name == "calcification rate" & standard_unit == "kg m^-2 year^-1", value/10, value),
         # standard_unit = ifelse(trait_name == "calcification rate" & standard_unit == "kg m^-2 year^-1", "g cm^-2 yr^-1", standard_unit),
         # # ^ changing value and standard unit for kg/m^2/yr -> g/cm^2/yr 
         # value = ifelse(trait_name == "calcification rate" & standard_unit == "mg CaCO3 cm^-2 d^-1", value*0.365, value),
         # standard_unit = ifelse(trait_name == "calcification rate" & standard_unit == "mg CaCO3 cm^-2 d^-1", "g cm^-2 yr^-1", standard_unit),
         # # ^ changing value and standard unit for mg CaCO3/cm^2/day -> g/cm^2/yr 
         # standard_unit = ifelse(trait_name == "calcification rate" & standard_unit == "g CaCO3 cm^-2 yr^-1", "g cm^-2 yr^-1", standard_unit),
         # # ^ changing unit only  (value stay the same)
         # standard_unit = ifelse(trait_name == "calcification rate" & standard_unit == "gCaCO3 cm^-2 yr^-1", "g cm^-2 yr^-1", standard_unit),
         # # ^ changing unit only  (value stay the same)
         # value = ifelse(trait_name == "calcification rate" & standard_unit == "kg CaCO3 m^-2", value/10, value), # I looked at primary literature, correct measurement is kg CaCO3/m^2/yr, year was not recorded on excel sheet (oops)
         # standard_unit = ifelse(trait_name == "calcification rate" & standard_unit == "kg CaCO3 m^-2", "g cm^-2 yr^-1", standard_unit),
         # # ^ changing value and standard unit for kg CaCO3/m^2/yr -> g/cm^2/yr 
         # value = ifelse(trait_name == "calcification rate" & standard_unit == "mmol CaCO3 m^-2 h^-1", value*0.08767884, value), 
         # standard_unit = ifelse(trait_name == "calcification rate" & standard_unit == "mmol CaCO3 m^-2 h^-1", "g cm^-2 yr^-1", standard_unit),
         # # ^ changing value and standard unit for mmol CaCO3/m^2/hr -> g/cm^2/yr
         # standard_unit = ifelse(trait_name == "calcification rate" & standard_unit == "g cm ^-2 yr^-1", "g cm^-2 yr^-1", standard_unit),
         # # ^ changing unit only  (value stay the same)
         # value = ifelse(trait_name == "calcification rate" & standard_unit == "mg CaCO3/cm^2/day", value*0.365, value),
         # standard_unit = ifelse(trait_name == "calcification rate" & standard_unit == "mg CaCO3/cm^2/day", "g cm^-2 yr^-1", standard_unit),
         # # ^ changing value and standard unit for mg CaCO3/cm^2/day -> g/cm^2/yr) 
         value = ifelse(trait_name == "symbiodinium density"  & standard_unit == "cells x 10^6 cm^-2", value/1e6, value),
         standard_unit = ifelse(trait_name == "symbiodinium density"  & standard_unit == "cells x 10^6 cm^-2", "cm^-2"   , standard_unit),
         # ^ changing value and standard unit for  "cells x 10^6 cm^-2" -> "cm^-2"
         # value = ifelse(trait_name == "chlorophyll a"  & standard_unit == "mg m^-2", value*0.1, value),
         # standard_unit = ifelse(trait_name == "chlorophyll a"  & standard_unit == "mg m^-2", "µg cm^-2"  , standard_unit),
         # # ^ changing value and standard unit for  "mg m^-2" -> "µg cm^-2"
  ) |> 
  filter(!(standard_unit %in% c("cm^2/month","g CaCO3 cm^-2", "% wt. gain month-1","Âµm cm^-2 h^-1"))) |> 
  select(specie_name, trait_name, value) |> 
  group_by(specie_name, trait_name) |> 
  summarize(mean = mean(value, na.rm = T),
            sd = sd(value, na.rm = T)) |> 
  group_by(specie_name) |> 
  mutate(sd2 = if_else(is.na(sd), mean(sd, na.rm = T), sd)) |> 
  ungroup() |> 
  select(-sd) |>
  pivot_wider(names_from = trait_name, values_from = c(mean, sd2))

  
# Log10-transform the trait values (excluding the species column)
log_traits <- as.data.frame(lapply(df_tr[,-1], function(x) log10(x)))

# Scale the log10-transformed trait values
scaled_traits <- as.data.frame(scale(log_traits))


# Combine the species column with the scaled trait values
scaled_traits <- cbind(df_tr[,1], scaled_traits)
colnames(scaled_traits)[1] <- colnames(df_tr)[1]

##load phylo tree 
phylo.tree<- funspace::phylo

coral_species <- phylo.tree$tip.label
print(class(coral_species))

missing_in_tree <- setdiff(coral_species, scaled_traits)
missing_in_data <- setdiff(scaled_traits, coral_species)

print("Species Missing in Tree:")
print(missing_in_tree)

print("Species Missing in Data:")
print(missing_in_data)

df_imputed <- impute(traits = scaled_traits, phylo = phylo.tree, addingSpecies = TRUE)



