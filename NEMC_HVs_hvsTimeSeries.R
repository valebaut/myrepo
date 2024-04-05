#' """ NEMC Time-Series Hypervolumes 
#'     @authors: Coral PR Team 
#'     date: """3/7/24

library(tidyverse)
library(hypervolume)
library(viridis)


# load trait data
df_tr = read_csv('zTraitTable.csv') |> 
  drop_na()

# load benthic data
df_ben = read.csv('benthicdata_PR.csv') %>% 
  filter(`SITE.NAME` %in% c("Carlos Rosario",  "Dakiti", "Cayo Diablo (2016)",  "Palominitos (2016)", "Palominos (2016)", "Canal Luis Pena" )) |> 
  filter(`YEAR` %in% c(2016, 2018, 2021, 2023)) %>% 
  group_by(species, `SITE.NAME`, YEAR) %>% 
  summarize(pc = mean(percentcover),
            pc_sd = sd(percentcover)) %>% 
  filter(pc>0)

#######make HVs 


reps = 100

set.seed(14)
df_NEMC = df_ben |> 
  left_join(df_tr, by = 'species') |> 
  ungroup() %>% 
  slice(rep(1:n(), each=reps))|> 
  mutate(i = rep(1:reps, times=nrow(df_ben)),
         percentcover = truncnorm::rtruncnorm(1, a = 0.001, b = 100,
                                              mean = pc, sd = pc_sd),
         `chlorophyll a` = truncnorm::rtruncnorm(1,
                                                 mean = `mean_chlorophyll a`,
                                                 sd = `sd2_chlorophyll a`),
         `corallite diameter` = truncnorm::rtruncnorm(1, 
                                                      mean = `mean_corallite diameter`,
                                                      sd = `sd2_corallite diameter`),
         `growth rate` = truncnorm::rtruncnorm(1, 
                                               mean = `mean_growth rate`,
                                               sd = `sd2_growth rate`),
         `skeletal density` = truncnorm::rtruncnorm(1, 
                                                    mean = `mean_skeletal density`,
                                                    sd = `sd2_skeletal density`),
         `symbiodinium density` = truncnorm::rtruncnorm(1, 
                                                        mean = `mean_symbiodinium density`,
                                                        sd = `sd2_symbiodinium density`),
         `colony maximum diameter` = truncnorm::rtruncnorm(1, 
                                                           mean = `mean_colony maximum diameter`,
                                                           sd = `sd2_colony maximum diameter`)) |> 
  select(`SITE.NAME`,YEAR, i:`colony maximum diameter`) |> 
  group_by(`SITE.NAME`,YEAR, i) |> 
  nest(weight = percentcover, data = `chlorophyll a`:`colony maximum diameter`) |> 
  mutate(hv = map2(data,weight, \(data,weight) hypervolume_gaussian(data, 
                                                                    name = paste(`SITE.NAME`, YEAR,i,sep = '_'),
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

saveRDS(df, 'NEMC_TimeSeries_hvs_randCov_randTr.rds')


df_NEMC$YEAR <- as.factor(df_NEMC$YEAR)
df_NEMC$'SITE.NAME' <- as.factor(df_NEMC$'SITE.NAME')



ggplot(df_NEMC, aes(YEAR, hv_size, color = SITE.NAME))+
  geom_point(size = 2.5)+
  geom_line(aes(group=SITE.NAME), linewidth = 1)+
  labs(x = 'Year', y = 'Volume')+
  theme_bw()+
  facet_wrap(~SITE.NAME, scales = 'free_y', nrow = 2)+
  #scale_y_log10()+
  #scale_color_viridis_d()+
  theme(axis.title = element_text(size = 14), 
        axis.text = element_text(size = 14, colour = "gray0"), 
        plot.title = element_text(size = 14, hjust=0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = 'right',
        legend.title = element_text(size = 14),
        strip.text.x = element_text(size = 14),
        legend.text = element_text(size = 12))

