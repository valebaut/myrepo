#' """ NEMC Time-Series BEM
#'     @authors: Coral PR Team 
#'     date: """4/5/24

library(tidyverse)
library(hypervolume)
library(viridis)


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
  select(YEAR, i:`colony maximum diameter`) |> 
  group_by(i) |> 
  mutate(across(`corallite diameter`:`colony maximum diameter`, scale)) |> 
  group_by(YEAR, i) |> 
  nest(weight = percentcover, data = `corallite diameter`:`colony maximum diameter`) |> 
  mutate(hv = map2(data,weight, \(data,weight) hypervolume_gaussian(data, 
                                                                    name = paste(YEAR,i,sep = '_'),
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


saveRDS(df, 'G:/R_analysis/YEAR_hvs_randCov_avgTr.rds')

###saving all the Hypervolumes 

df |> select(YEAR, i, hv_size, centroid) |> 
  unnest_wider(centroid) |> 
  write_csv('G:/R_analysis/YEAR_hvs_randCov_avgTr.csv')
# will be too big for github

# # plot size across YEAR
# ggplot(df, aes(YEAR, hv_size))+
#   geom_boxplot()


# comparison of across YEARS (change df_y to gb)
df_reg = tibble(YEAR1 = unique(df$YEAR),
                YEAR2 = unique(df$YEAR)) |> 
  expand(YEAR1, YEAR2)

df_reg = df_reg[!duplicated(t(apply(df_reg,1,sort))),] %>% 
  filter(!(YEAR1 == YEAR2))

df_iterations = tibble(YEAR1 = rep(df_reg$YEAR1, times = reps),
                       YEAR2 = rep(df_reg$YEAR2, times = reps),
                       i = rep(1:reps, each = length(df_reg$YEAR1))) 

# years to make 
df_1 = df|> 
  select(YEAR1 = YEAR, hv1 = hv, hv1_size = hv_size, i = i)

df_2 = df |> 
  select(YEAR2 = YEAR, hv2 = hv, hv2_size = hv_size, i = i)

# create large df to store all data

# tibble(BASIN = rep(unique(df$BASIN),
#             each = nrow(df_y)),
# y1 = rep(df_y$y1, times = length(unique(df$BASIN))),
# y2 = rep(df_y$y2, times = length(unique(df$BASIN))))

df_ov = df_iterations |> 
  inner_join(df_1, by = c('YEAR1', 'i')) |> 
  inner_join(df_2, by = c('YEAR2', 'i')) |> 
  mutate(set = map2(hv1,hv2, \(hv1, hv2) hypervolume_set(hv1, hv2, check.memory = F, verbose = F)),
         ov = map(set, \(set) hypervolume_overlap_statistics(set)),
         dist_cent = map2_dbl(hv1, hv2, \(hv1,hv2) hypervolume_distance(hv1, hv2, type = 'centroid', check.memory=F)),
         size_ratio = hv1_size/hv2_size,
         i = i) |> 
  unnest_wider(ov) |> 
  select(YEAR1, YEAR2, hv1_size, hv2_size, size_ratio,
         jaccard, sorensen,uniq_y1 = frac_unique_1, uniq_y2 = frac_unique_2, 
         dist_cent)

#saveRDS(df_ov, "data/YEAR_ov.rds")
write_csv(df_ov, "G:/R_analysis/YEAR_ov_randCov_avgTr.csv")

df_avg = df_ov |> 
  group_by(YEAR1,YEAR2) |>
  summarize(across(hv1_size:dist_cent, mean))

ggplot(df_avg, aes(YEAR1, YEAR2, fill = sorensen))+
  geom_tile()+
  scale_fill_viridis(limits = c(0,1))

