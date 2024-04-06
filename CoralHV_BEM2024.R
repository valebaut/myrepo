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
reps = 100

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


saveRDS(df, 'C:/Users/valeb/OneDrive - Florida International University/GitHubRep/NEMC_hvs_randCov_avgTr_.rds')

###saving all the Hypervolumes 

df |> 
  select(YEAR, site,i, hv_size, centroid) |> 
  unnest_wider(centroid) |> 
  write_csv('hvALL.csv')

# will be too big for github

# # plot size across YEAR
# ggplot(df, aes(YEAR, hv_size))+
#   geom_boxplot()

##############SITES OVER TIME 
df_y= tibble(y1 = unique(df$YEAR),
             y2 = unique(df$YEAR)) |> 
  expand(y1,y2)

df_y = df_y[!duplicated(t(apply(df_y,1,sort))),] %>% 
  filter(!(y1 == y2))

# years to make 
df1 = df |> 
  select(site, y1 = YEAR, hv1 = hv, hv1_size = hv_size, cent1 = centroid, i)

df2 = df |> 
  select(site, y2 = YEAR, hv2 = hv, hv2_size = hv_size, cent2 = centroid, i)


# create data frame of all data and make yearly comparisons
sitesthrutime = tibble(site = rep(unique(df$site),
                           each = nrow(df_y)),
               y1 = rep(df_y$y1, times = length(unique(df$site))),
               y2 = rep(df_y$y2, times = length(unique(df$site)))) |> 
  inner_join(df1, by = c('site', 'y1')) |> 
  inner_join(df2, by = c('site', 'y2')) |> 
  mutate(ychange = y2-y1,
         size_rat = hv2_size/hv1_size,
         set = map2(hv1,hv2, \(hv1, hv2) hypervolume_set(hv1, hv2, check.memory = F, verbose = F)),
         ov = map(set, \(set) hypervolume_overlap_statistics(set)),
         dist_cent = map2_dbl(hv1, hv2, \(hv1,hv2) hypervolume_distance(hv1, hv2, type = 'centroid', check.memory=F)),
         dif = map2(cent1, cent2, \(cent1,cent2) cent2 - cent1)) |> 
  unnest_wider(ov) |> 
  unnest_wider(dif) |> 
  select(site, y1, y2, ychange, hv1_size, hv2_size, size_rat, 
         jaccard, sorensen,uniq_y1 = frac_unique_1, uniq_y2 = frac_unique_2, 
         dist_cent)

write_csv(sitesthrutime, "sitesthrutime.csv")

##########WITHIN YEAR COMPARISON ACROSS SITES 
####seperate years 
WI_2016 = df %>% 
  filter(YEAR==2016)
WI_2018 = df %>% 
  filter(YEAR==2018)
WI_2021 = df %>% 
  filter(YEAR==2021)
WI_2023 = df %>% 
  filter(YEAR==2023)
####build tibble for each year 
##2016
reg_2016 = tibble(site1 = unique(WI_2016$site),
                site2 = unique(WI_2016$site)) |> 
  expand(site1,site2)

reg_2016 = reg_2016[!duplicated(t(apply(reg_2016,1,sort))),] %>% 
  filter(!(site1 == site2))

iteration_2016 = tibble(site1 = rep(reg_2016$site1, times = reps),
                       site2 = rep(reg_2016$site2, times = reps),
                       i = rep(1:reps, each = length(reg_2016$site1))) 

# years to make 
df_161 = WI_2016|> 
  select(site1 = site, hv1 = hv, hv1_size = hv_size, i = i, YEAR)

df_162 = WI_2016 |> 
  select(site2 = site, hv2 = hv, hv2_size = hv_size, i = i, YEAR)
###2018
reg_2018 = tibble(site1 = unique(WI_2018$site),
                  site2 = unique(WI_2018$site)) |> 
  expand(site1,site2)

reg_2018 = reg_2018[!duplicated(t(apply(reg_2018,1,sort))),] %>% 
  filter(!(site1 == site2))

iteration_2018 = tibble(site1 = rep(reg_2018$site1, times = reps),
                        site2 = rep(reg_2018$site2, times = reps),
                        i = rep(1:reps, each = length(reg_2018$site1))) 

# years to make 
df_181 = WI_2018|> 
  select(site1 = site, hv1 = hv, hv1_size = hv_size, i = i, YEAR)

df_182 = WI_2018 |> 
  select(site2 = site, hv2 = hv, hv2_size = hv_size, i = i, YEAR)
###2021
reg_2021 = tibble(site1 = unique(WI_2021$site),
                  site2 = unique(WI_2021$site)) |> 
  expand(site1,site2)

reg_2021 = reg_2021[!duplicated(t(apply(reg_2021,1,sort))),] %>% 
  filter(!(site1 == site2))

iteration_2021 = tibble(site1 = rep(reg_2021$site1, times = reps),
                        site2 = rep(reg_2021$site2, times = reps),
                        i = rep(1:reps, each = length(reg_2021$site1))) 

# years to make 
df_211 = WI_2021|> 
  select(site1 = site, hv1 = hv, hv1_size = hv_size, i = i, YEAR)

df_212 = WI_2021 |> 
  select(site2 = site, hv2 = hv, hv2_size = hv_size, i = i, YEAR)
###2023
reg_2023 = tibble(site1 = unique(WI_2023$site),
                  site2 = unique(WI_2023$site)) |> 
  expand(site1,site2)

reg_2023 = reg_2023[!duplicated(t(apply(reg_2023,1,sort))),] %>% 
  filter(!(site1 == site2))

iteration_2023 = tibble(site1 = rep(reg_2023$site1, times = reps),
                        site2 = rep(reg_2023$site2, times = reps),
                        i = rep(1:reps, each = length(reg_2023$site1))) 

# years to make 
df_231 = WI_2023|> 
  select(site1 = site, hv1 = hv, hv1_size = hv_size, i = i, YEAR)

df_232 = WI_2023 |> 
  select(site2 = site, hv2 = hv, hv2_size = hv_size, i = i, YEAR)
####Make iterations for each year 
###2016
ov_2016 = iteration_2016 |> 
  inner_join(df_161, by = c('site1', 'i')) |> 
  inner_join(df_162, by = c('site2', 'i')) |> 
  mutate(set = map2(hv1,hv2, \(hv1, hv2) hypervolume_set(hv1, hv2, check.memory = F, verbose = F)),
         ov = map(set, \(set) hypervolume_overlap_statistics(set)),
         dist_cent = map2_dbl(hv1, hv2, \(hv1,hv2) hypervolume_distance(hv1, hv2, type = 'centroid', check.memory=F)),
         size_ratio = hv1_size/hv2_size,
         i = i) |> 
  mutate(year=2016) |>
  unnest_wider(ov) |> 
  select(site1, site2, hv1_size, hv2_size, size_ratio,
         jaccard, sorensen,uniq_y1 = frac_unique_1, uniq_y2 = frac_unique_2, 
         dist_cent,i, year)

###2018 
ov_2018 = iteration_2018 |> 
  inner_join(df_181, by = c('site1', 'i')) |> 
  inner_join(df_182, by = c('site2', 'i')) |> 
  mutate(set = map2(hv1,hv2, \(hv1, hv2) hypervolume_set(hv1, hv2, check.memory = F, verbose = F)),
         ov = map(set, \(set) hypervolume_overlap_statistics(set)),
         dist_cent = map2_dbl(hv1, hv2, \(hv1,hv2) hypervolume_distance(hv1, hv2, type = 'centroid', check.memory=F)),
         size_ratio = hv1_size/hv2_size,
         i = i) |> 
  mutate(year=2018) |>
  unnest_wider(ov) |> 
  select(site1, site2, hv1_size, hv2_size, size_ratio,
         jaccard, sorensen,uniq_y1 = frac_unique_1, uniq_y2 = frac_unique_2, 
         dist_cent,i, year)
###2021 
ov_2021 = iteration_2021 |> 
  inner_join(df_211, by = c('site1', 'i')) |> 
  inner_join(df_212, by = c('site2', 'i')) |> 
  mutate(set = map2(hv1,hv2, \(hv1, hv2) hypervolume_set(hv1, hv2, check.memory = F, verbose = F)),
         ov = map(set, \(set) hypervolume_overlap_statistics(set)),
         dist_cent = map2_dbl(hv1, hv2, \(hv1,hv2) hypervolume_distance(hv1, hv2, type = 'centroid', check.memory=F)),
         size_ratio = hv1_size/hv2_size,
         i = i) |> 
  mutate(year=2021) |>
  unnest_wider(ov) |> 
  select(site1, site2, hv1_size, hv2_size, size_ratio,
         jaccard, sorensen,uniq_y1 = frac_unique_1, uniq_y2 = frac_unique_2, 
         dist_cent,i, year)
###2022 
ov_2023 = iteration_2023 |> 
  inner_join(df_231, by = c('site1', 'i')) |> 
  inner_join(df_232, by = c('site2', 'i')) |> 
  mutate(set = map2(hv1,hv2, \(hv1, hv2) hypervolume_set(hv1, hv2, check.memory = F, verbose = F)),
         ov = map(set, \(set) hypervolume_overlap_statistics(set)),
         dist_cent = map2_dbl(hv1, hv2, \(hv1,hv2) hypervolume_distance(hv1, hv2, type = 'centroid', check.memory=F)),
         size_ratio = hv1_size/hv2_size,
         i = i) |> 
  mutate(year=2022) |>
  unnest_wider(ov) |> 
  select(site1, site2, hv1_size, hv2_size, size_ratio,
         jaccard, sorensen,uniq_y1 = frac_unique_1, uniq_y2 = frac_unique_2, 
         dist_cent,i, year)

###Join all four years 
joined <- bind_rows(ov_2016, ov_2018, ov_2021, ov_2023) 

#saveRDS(df_ov, "data/site_ov.rds")
write_csv(joined, "withinyear_acrosssites.csv")

################################################################



df_avg = df_ov |> 
  group_by(YEAR1,YEAR2) |>
  summarize(across(hv1_size:dist_cent, mean))

ggplot(df_avg, aes(YEAR1, YEAR2, fill = sorensen))+
  geom_tile()+
  scale_fill_viridis(limits = c(0,1))

