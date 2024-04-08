#' """ Quantify responses
#'     @author: Ryan James
#'     Date : 10/12/23

library(tidyverse)
library(performance)
library(MuMIn)
library(hypervolume)
library(viridis)

# step 1 variance ----
##Within site change over time 
# centroid distance per site through time 
site_time = read_csv("sitesthrutime.csv")  
df_var <- site_time %>%
  group_by(site, y1, y2) %>%
  reframe(
    var = var(dist_cent),
    mean = mean(dist_cent),
    sd = sd(dist_cent),
    cv = sd / mean,
    stab = 1 / cv,
    ychange = unique(ychange)
  )
unique(df_var$site)
df_var$site <- factor(df_var$site, levels = c("Canal Luis Pena", "Carlos Rosario", "Dakiti", "Cayo Diablo (2016)", "Palominitos (2016)", "Palominos (2016)"))

ggplot(df_var, aes(ychange, mean, color = site))+
    geom_hline(aes(yintercept = 1), linetype = 'dashed')+
    geom_point(size = 2.5)+
    geom_smooth(method = 'lm', formula = y~x, 
              linewidth = 1, color = 'black')+
    facet_wrap(~site,  nrow = 2)+
    scale_color_viridis_d(option = 'turbo')+
    labs(x = 'Years Between', y = 'Centroid distance')+theme_bw()+
    theme(axis.title = element_text(size = 14), 
        axis.text.y = element_text(size = 14, colour = "black"),
        axis.text.x = element_text(size = 12, colour = "black"), 
        plot.title = element_text(size = 14, hjust=0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = 'none',
        legend.title = element_text(size = 14),
        strip.text.x = element_text(size = 14),
        legend.text = element_text(size = 12)) 
ggsave('figs/centResp.png', 
       units="in", width=10, height=6, dpi=600)


# overlap per site through time 
site_time = read_csv("sitesthrutime.csv")  
df_var <- site_time %>%
  group_by(site, y1, y2) %>%
  reframe(
    var = var(sorensen),
    mean = mean(sorensen),
    sd = sd(sorensen),
    cv = sd / mean,
    stab = 1 / cv,
    ychange = unique(ychange)
  )

df_var$site <- factor(df_var$site, levels = c("Canal Luis Pena", "Carlos Rosario", "Dakiti", "Cayo Diablo (2016)", "Palominitos (2016)", "Palominos (2016)"))
ggplot(df_var, aes(ychange, mean, color = site))+
  geom_hline(aes(yintercept = 0), linetype = 'dashed')+
  geom_point(size = 2.5)+
  geom_smooth(method = 'lm', formula = y~x, 
              linewidth = 1, color = 'black')+
  facet_wrap(~site,  nrow = 2)+
  scale_color_viridis_d(option = 'turbo')+
  labs(x = 'Years Between', y = 'Overlap')+theme_bw()+
  theme(axis.title = element_text(size = 14), 
        axis.text.y = element_text(size = 14, colour = "black"),
        axis.text.x = element_text(size = 12, colour = "black"), 
        plot.title = element_text(size = 14, hjust=0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = 'none',
        legend.title = element_text(size = 14),
        strip.text.x = element_text(size = 14),
        legend.text = element_text(size = 12)) 
ggsave('figs/overlapResp.png', 
       units="in", width=10, height=6, dpi=600)

# volume per site through time 
df_var <- site_time %>%
  group_by(site, y1, y2) %>%
  reframe(
    var = var(hv1_size),
    mean = mean(hv1_size),
    sd = sd(hv1_size),
    cv = sd / mean,
    stab = 1 / cv,
    ychange = unique(ychange)
  )
df_var$site <- factor(df_var$site, levels = c("Canal Luis Pena", "Carlos Rosario", "Dakiti", "Cayo Diablo (2016)", "Palominitos (2016)", "Palominos (2016)"))
ggplot(df_var, aes(y1, mean, color = site))+
  geom_hline(aes(yintercept = 0), linetype = 'dashed')+
  geom_point(size = 2.5)+
  geom_smooth(method = 'lm', formula = y~x, 
              linewidth = 1)+
  facet_wrap(~site,  nrow = 2)+
  labs(x = 'Year', y = 'Volume')+
  scale_color_viridis_d(option = 'turbo')+
  #scale_y_log10()+
  theme_bw()+
  theme(axis.title = element_text(size = 14),
        axis.text.y = element_text(size = 14, colour = "black"),
        axis.text.x = element_text(size = 12, colour = "black"), 
        plot.title = element_text(size = 14, hjust=0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = 'none',
        legend.title = element_text(size = 14),
        strip.text.x = element_text(size = 14),
        legend.text = element_text(size = 12))

ggsave('figs/sizeResp.png', 
       units="in", width=10, height=6, dpi=600)

# size ratio----
df = read_csv('hvALL.csv')

df_ov = read_csv("sitesthrutime.csv") |> 
  mutate(lsr = log(size_rat)) |>
  group_by(site, y1, y2) |> 
  summarise(
    lsr = mean(lsr),
    sd = sd(lsr),
    ychange = unique(ychange))|> 
  group_by(site) |>
  nest() |> 
  mutate(m_int = map(data, \(df)lm(lsr~1, data = df)),
         m_lin = map(data, \(df)lm(lsr~ychange, data = df)),
         m_quad = map(data, \(df)lm(lsr~ychange + I(ychange^2), data = df)),
         AICc_int = map_dbl(m_int, \(x) AICc(x)),
         AICc_lin = map_dbl(m_lin, \(x) AICc(x)),
         AICc_quad = map_dbl(m_quad, \(x) AICc(x)),
         model = case_when(
           AICc_int - min(c(AICc_int,AICc_lin,AICc_quad)) <= 4 ~ 'Intercept',
           AICc_lin < AICc_quad ~ 'Linear',
           AICc_quad < AICc_lin ~ 'Quadratic'))

d = df_ov |> 
  select(site, data, model) |> 
  unnest(cols = c(data))
d$site <- factor(d$site, levels = c("Canal Luis Pena", "Carlos Rosario", "Dakiti", "Cayo Diablo (2016)", "Palominitos (2016)", "Palominos (2016)"))

ggplot(d, aes(ychange, lsr, color = site))+
  geom_hline(aes(yintercept = 0), linetype = 'dashed')+
  geom_point(size = 2.5)+
  # geom_smooth(data = d |> filter(model == 'Intercept'),
  #             method = 'lm', formula = y~1, 
  #             linewidth = 1, color = 'black')+
  geom_smooth(method = 'lm', formula = y~x, 
              linewidth = 1, color = 'black')+
  # geom_smooth(data = d |> filter(model == 'Quadratic'),
  #             method = 'lm', formula = y~x+I(x^2), 
  #             linewidth = 1, color = 'black')+
  facet_wrap(~site,  nrow = 2)+
  scale_color_viridis_d(option = 'turbo')+
  labs(x = 'Years between comparison', y = 'log(Y2/Y1 size ratio)')+theme_bw()+
  theme(axis.title = element_text(size = 14), 
        axis.text.y = element_text(size = 14, colour = "black"),
        axis.text.x = element_text(size = 12, colour = "black"), 
        plot.title = element_text(size = 14, hjust=0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = 'none',
        legend.title = element_text(size = 14),
        strip.text.x = element_text(size = 14),
        legend.text = element_text(size = 12))

ggsave('figs/sizeRespLog.png', 
       units="in", width=10, height=6, dpi=600)



###Across sites within a single year  
###centroid distance 
withiyear_acrosssite = read_csv("withinyear_acrosssites.csv")  %>%
  group_by(site1, site2, year) %>%
  mutate(year = if_else(year == 2022, 2023, year)) |> 
  reframe(
    var = var(dist_cent),
    mean = mean(dist_cent),
    sd = sd(dist_cent),
    cv = sd / mean,
    stab = 1 / cv
  ) 
# |> 
#   group_by(year) |> 
#   summarise(
#     mean=mean(mean)
#   )


ggplot(withiyear_acrosssite, aes(year, mean))+
  geom_hline(aes(yintercept = 0), linetype = 'dashed')+
  geom_point(size = 2.5)+
  geom_smooth(method = 'lm', formula = y~x, 
              linewidth = 1, color = 'black')+
  #facet_wrap(~site1,  nrow = 2)+
  scale_color_viridis_d(option = 'turbo')+
  labs(x = 'Year', y = 'Centroid distance')+theme_bw()+
  theme(axis.title = element_text(size = 14), 
        axis.text.y = element_text(size = 14, colour = "black"),
        axis.text.x = element_text(size = 12, colour = "black"), 
        plot.title = element_text(size = 14, hjust=0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = 'none',
        legend.title = element_text(size = 14),
        strip.text.x = element_text(size = 14),
        legend.text = element_text(size = 12)) 
ggsave('figs/hvDistYearly.png', 
       units="in", width=10, height=6, dpi=600)
##overlap
withiyear_acrosssite = read_csv("withinyear_acrosssites.csv")  %>%
  group_by(site1, site2, year) %>%
  mutate(year = if_else(year == 2022, 2023, year)) |> 
  reframe(
    var = var(sorensen),
    mean = mean(sorensen),
    sd = sd(sorensen),
    cv = sd / mean,
    stab = 1 / cv
   ) #|> 
  # group_by(year) |> 
  # summarise(
  #   mean=mean(mean)
  # )


ggplot(withiyear_acrosssite, aes(year, mean))+
  geom_hline(aes(yintercept = 0), linetype = 'dashed')+
  geom_point(size = 2.5)+
  geom_smooth(method = 'lm', formula = y~x, 
              linewidth = 1, color = 'black')+
  #facet_wrap(~site1,  nrow = 2)+
  scale_color_viridis_d(option = 'turbo')+
  labs(x = 'Year', y = 'Overlap')+theme_bw()+
  theme(axis.title = element_text(size = 14), 
        axis.text.y = element_text(size = 14, colour = "black"),
        axis.text.x = element_text(size = 12, colour = "black"), 
        plot.title = element_text(size = 14, hjust=0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = 'none',
        legend.title = element_text(size = 14),
        strip.text.x = element_text(size = 14),
        legend.text = element_text(size = 12)) 

ggsave('figs/hvoverlapYearly.png', 
       units="in", width=10, height=6, dpi=600)
###size/ volume 
withiyear_acrosssite = read_csv('hvALL.csv')  %>%
  group_by(site, YEAR) %>%
  reframe(
    var = var(hv_size),
    mean = mean(hv_size),
    sd = sd(hv_size),
    cv = sd / mean,
    stab = 1 / cv
   )  #|> 
  # group_by(YEAR) |> 
  # summarise(
  #   mean=mean(mean)
  # )

ggplot(withiyear_acrosssite, aes(YEAR, mean))+
  geom_hline(aes(yintercept = 0), linetype = 'dashed')+
  geom_point(size = 2.5)+
  geom_smooth(method = 'lm', formula = y~x, 
              linewidth = 1, color = 'black')+
  #facet_wrap(~site1,  nrow = 2)+
  scale_color_viridis_d(option = 'turbo')+
  labs(x = 'Year', y = 'Trait Space Size')+theme_bw()+
  theme(axis.title = element_text(size = 14), 
        axis.text.y = element_text(size = 14, colour = "black"),
        axis.text.x = element_text(size = 12, colour = "black"), 
        plot.title = element_text(size = 14, hjust=0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = 'none',
        legend.title = element_text(size = 14),
        strip.text.x = element_text(size = 14),
        legend.text = element_text(size = 12)) 
ggsave('figs/hvsizeYearly.png', 
       units="in", width=10, height=6, dpi=600)

##Step 2: trait importance
#####trait centroid 
withiyear_acrosssite = read_csv('hvALL.csv')  |> 
  pivot_longer(corallite.diameter:colony.maximum.diameter, names_to = "traits", values_to = "value") |> 
  group_by(site, YEAR, traits) |> 
  reframe(
    var = var(value),
    mean = mean(value),
    sd = sd(value),
    cv = sd / mean,
    stab = 1 / cv
  ) |> 
  group_by(YEAR, traits) |> 
  summarise(
    mean=mean(mean)
  )

ggplot(withiyear_acrosssite, aes(YEAR, mean, color=traits))+
  geom_hline(aes(yintercept = 0), linetype = 'dashed')+
  geom_point(size = 2.5)+
  geom_smooth(method = 'lm', formula = y~x, 
              linewidth = 1)+
  facet_wrap(~traits,  nrow = 2)+
  scale_color_viridis_d(option = 'turbo')+
  labs(x = 'Year', y = 'Trait Centroid')+theme_bw()+
  theme(axis.title = element_text(size = 14), 
        axis.text.y = element_text(size = 14, colour = "black"),
        axis.text.x = element_text(size = 12, colour = "black"), 
        plot.title = element_text(size = 14, hjust=0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = 'none',
        legend.title = element_text(size = 14),
        strip.text.x = element_text(size = 14),
        legend.text = element_text(size = 12)) 

ggsave('figs/TraitCentroidYearly.png', 
       units="in", width=10, height=6, dpi=600)



