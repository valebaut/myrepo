#' """ NEMC Time-Series BEM
#'     @authors: Coral PR Team 
#'     date: """4/5/24

library(tidyverse)
library(hypervolume)
library(viridis)
library(purrr)
library(dplyr)

# load trait data
df_tr = read_csv('TraitTable.csv') |> 
  drop_na()

# load benthic data
df_ben = read.csv('benthicdata_PR.csv')|> 
rename(site= SITE.NAME) |>
  filter(`site` %in% c("Carlos Rosario",  "Dakiti", "Cayo Diablo (2016)",  "Palominitos (2016)", "Palominos (2016)", "Canal Luis Pena" )) |> 
  filter(`YEAR` %in% c(2016, 2018, 2021, 2023))|>
  filter(species %in% unique(df_tr$species)) |> 
  pivot_wider(names_from =species, values_from = percentcover, values_fill = 0) |> 
  pivot_longer(`Acropora cervicornis`:`Orbicella franksi`, 
               names_to = 'species', values_to = 'percentcover') |> 
  group_by(species, YEAR, site) |> 
  summarize(pc = mean(percentcover),
            pc_sd = sd(percentcover),
            .groups = 'drop') |> 
  filter(pc > 0) 

# make the hvs for random cover and avg traits----
reps = 1

set.seed(14)
df = df_ben |> 
  left_join(df_tr, by = 'species') |> 
  slice(rep(1:n(), each=reps))|> 
  mutate(i = rep(1:reps, times=nrow(df_ben)),
         percentcover = truncnorm::rtruncnorm(1, a = 0.001, b = 100,
                                              mean = pc, sd = pc_sd),
         # `chlorophyll a` = truncnorm::rtruncnorm(1,
         #                                              mean = `mean_chlorophyll a`,
         #                                              sd = `sd2_chlorophyll a`),
         `corallite diameter` = `mean_corallite diameter`,
         `growth rate` = `mean_growth rate`,
         `skeletal density` = `mean_skeletal density`,
         `symbiodinium density` = `mean_symbiodinium density`,
         `colony maximum diameter` = `mean_colony maximum diameter`) |> 
  select(YEAR, site, i:`colony maximum diameter`) |> 
  group_by(i) |> 
  mutate(across(`corallite diameter`:`colony maximum diameter`, scale)) |> 
  group_by(YEAR, i, site) |> 
  nest(weight = percentcover, data = `corallite diameter`:`colony maximum diameter`) |> 
  mutate(hv = map2(data,weight, \(data,weight) hypervolume_gaussian(data, 
                                                                    name = paste(YEAR,site,i,sep = '_'),
                                                                    weight = weight$percentcover,
                                                                    samples.per.point = 1000,
                                                                    kde.bandwidth = estimate_bandwidth(data), 
                                                                    sd.count = 3, 
                                                                    quantile.requested = 0.95, 
                                                                    quantile.requested.type = "probability", 
                                                                    chunk.size = 1000, 
                                                                    verbose = F)),
         hv_size = map_dbl(hv, \(hv) get_volume(hv)),
         centroid = map(hv, \(hv) get_centroid(hv)))


saveRDS(df, 'C:/Users/valeb/OneDrive - Florida International University/GitHubRep/YEAR_hvs_randCov_avgTr_SITES.rds')

###saving all the Hypervolumes 

df |> select(YEAR, site,i, hv_size, centroid) |> 
  unnest_wider(centroid) |> 
  write_csv('YEAR_hvs_randCov_avgTr_SITES.csv')
# will be too big for github

# # plot size across YEAR
# ggplot(df, aes(YEAR, hv_size))+
#   geom_boxplot()

##############SITES OVER TIME 

# Extract unique values of YEAR and site from the dataframe df
unique_years <- unique(df$YEAR)
unique_sites <- unique(df$site)

# Create tibble with combinations of YEAR1, YEAR2, site1, and site2
df_reg <- expand.grid(YEAR1 = unique_years, YEAR2 = unique_years,
                      site1 = unique_sites, site2 = unique_sites) %>%
  filter(YEAR1 != YEAR2)  # Filter out combinations where YEAR1 equals YEAR2


df_reg = df_reg[!duplicated(t(apply(df_reg,1,sort))),] 

df_iterations = tibble(YEAR1 = rep(df_reg$YEAR1, times = reps),
                       YEAR2 = rep(df_reg$YEAR2, times = reps),
                       site1 = rep(df_reg$site1, times = reps),  # Adding site1
                       site2 = rep(df_reg$site2, times = reps),  # Adding site2
                       i = rep(1:reps, each = length(df_reg$YEAR1))) 

df_1 = df |> 
  select(YEAR1 = YEAR, site1 = site, hv1 = hv, hv1_size = hv_size, i = i)  # Adding site1

df_2 = df |> 
  select(YEAR2 = YEAR, site2 = site, hv2 = hv, hv2_size = hv_size, i = i)  # Adding site2

df_ov = df_iterations |> 
  inner_join(df_1, by = c('YEAR1', 'site1', 'i')) |>  # Adding site1
  inner_join(df_2, by = c('YEAR2', 'site2', 'i')) |>  # Adding site2
  mutate(set = map2(hv1, hv2, \(hv1, hv2) hypervolume_set(hv1, hv2, check.memory = F, verbose = F)),
         ov = map(set, \(set) hypervolume_overlap_statistics(set)),
         dist_cent = map2_dbl(hv1, hv2, \(hv1, hv2) hypervolume_distance(hv1, hv2, type = 'centroid', check.memory = F)),
         size_ratio = hv1_size/hv2_size,
         i = i) |> 
  unnest_wider(ov) |> 
  select(YEAR1, YEAR2, site1, site2, hv1_size, hv2_size, size_ratio,
         jaccard, sorensen, uniq_y1 = frac_unique_1, uniq_y2 = frac_unique_2, 
         dist_cent)


#saveRDS(df_ov, "data/YEAR_ov.rds")
write_csv(df_ov, "YEAR_ov_randCov_avgTr_SITES.csv")



##########WITHIN YEAR COMPARISON ACROSS SITES 

# Generate combinations of sites within the same year
df_site_combinations <- df %>%
  group_by(YEAR) %>%
  mutate(combinations = list(combn(site, 2, simplify = FALSE))) %>%
  unnest(combinations) %>%
  rename(SITE1 = combinations[[1]], SITE2 = combinations[[2]]) %>%
  distinct()

# Perform comparisons between sites within the same year
df_ov_2 <- df_site_combinations %>%
  mutate(i = row_number()) %>%
  left_join(df, by = c("YEAR", "SITE1")) %>%
  select(-site) %>%
  rename(hv1 = hv, hv1_size = hv_size) %>%
  left_join(df, by = c("YEAR", "SITE2")) %>%
  select(-site) %>%
  rename(hv2 = hv, hv2_size = hv_size) %>%
  mutate(set = map2(hv1, hv2, \(hv1, hv2) hypervolume_set(hv1, hv2, check.memory = FALSE, verbose = FALSE)),
         ov = map(set, \(set) hypervolume_overlap_statistics(set)),
         dist_cent = map2_dbl(hv1, hv2, \(hv1, hv2) hypervolume_distance(hv1, hv2, type = 'centroid', check.memory = FALSE)),
         size_ratio = hv1_size / hv2_size) %>%
  unnest_wider(ov) %>%
  select(YEAR, site1, SITE2, hv1_size, hv2_size, size_ratio, jaccard, sorensen,
         uniq_y1 = frac_unique_1, uniq_y2 = frac_unique_2, dist_cent)

# Save the resulting dataframe
write.csv(df_ov_2, "comparison_within_year.csv", row.names = FALSE)





withinyear<-df_1 %>% 
  select(YEAR1, site1, hv_size) |> 
  unnest_wider(centroid) |> 
  write_csv(withinyear, "YEAR_ov_randCov_avgTr_withinyear.csv")

df_avg = df_ov |> 
  group_by(YEAR1,YEAR2) |>
  summarize(across(hv1_size:dist_cent, mean))

ggplot(df_avg, aes(YEAR1, YEAR2, fill = sorensen))+
  geom_tile()+
  scale_fill_viridis(limits = c(0,1))

