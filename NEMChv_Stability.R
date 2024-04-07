#' """ Quantify responses
#'     @author: Ryan James
#'     Date : 10/12/23

library(tidyverse)
library(performance)
library(MuMIn)
library(hypervolume)
library(viridis)

# step 1 variance ----
# plot centroid dist and size over year ----
# centroid distance
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


###site over time 
ggplot(df_var, aes(ychange, mean, color = site))+
    geom_hline(aes(yintercept = 1), linetype = 'dashed')+
    geom_point(size = 2.5)+
    geom_line(linewidth = 1)+
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

###within year comparison across sites 
withiyear_acrosssite = read_csv("withinyear_acrosssites.csv")  %>%
  group_by(site1, site2, year) %>%
  reframe(
    var = var(dist_cent),
    mean = mean(dist_cent),
    sd = sd(dist_cent),
    cv = sd / mean,
    stab = 1 / cv
  )

ggplot(withiyear_acrosssite, aes(year, mean, color = site1))+
  geom_hline(aes(yintercept = 1), linetype = 'dashed')+
  geom_point(size = 2.5)+
  geom_line(linewidth = 1)+
  facet_wrap(~site1,  nrow = 2)+
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

# size 
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

ggplot(df_var, aes(y1, mean, color = site))+
  geom_point(size = 2.5)+
  geom_line(linewidth = 1)+
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

ggsave('figs/hvSizeYearly.png', 
       units="in", width=10, height=6, dpi=600)

# centroid distance loo ----
df_ov = read_csv("sitesthrutime.csv") |> 
  filter(ychange == 2)

ys = unique(df_ov[c('y1','y2')])

for(i in 1:nrow(ys)){
  d = df_ov |> 
    #filter(y1 != ys$y1[i], y2 != ys$y2[i]) |> 
    group_by(site) |> 
    summarize(var = var(dist_cent),
              mean = mean(dist_cent),
              sd = sd(dist_cent),
              cv = sd/mean,
                stab = 1/cv) 
  
  if(i == 1){
    df_var = d
  }else{
    df_var = bind_rows(df_var, d)
  }}


df_mm = df_var |> 
  group_by(site) |> 
  summarize(lc = quantile(stab, 0.025),
            uc = quantile(stab, 0.975))

df_cdv = df_ov |> 
  group_by(site) |> 
  summarize(var = var(dist_cent),
            mean = mean(dist_cent),
            sd = sd(dist_cent),
            cv = sd/mean,
            stab = 1/cv) |> 
  left_join(df_mm, by = 'site') 

ggplot(df_cdv, aes(site, stab, color = site))+
  geom_point(size = 5)+
  geom_errorbar(aes(ymin = lc, ymax = uc), linewidth = 2, width = 0)+
  labs(x = 'site', y = 'Stability centroid distance')+
  scale_color_viridis_d(option = 'turbo')+
  theme_bw()+
  theme(axis.title = element_text(size = 14), 
        axis.text = element_text(size = 14, colour = "gray0"), 
        plot.title = element_text(size = 14, hjust=0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = 'none',
        legend.title = element_text(size = 14),
        strip.text.x = element_text(size = 14),
        legend.text = element_text(size = 12))

ggsave('figs/stabCentDist.png', 
       units="in", width=8, height=6, dpi=600)


# hv size loo ----
df = read_csv('hvALL.csv') 

ys = unique(df$YEAR)

for(i in 1:length(ys)){
  d = df |> 
    filter(YEAR != ys[i]) |> 
    group_by(site) |> 
    summarize(var = var(hv_size),
              mean = mean(hv_size),
              sd = sd(hv_size),
              cv = sd/mean,
              stab = 1/cv)
  
  if(i == 1){
    df_var = d
  }else{
    df_var = bind_rows(df_var, d)
  }
}

df_mm = df_var |> 
  group_by(site) |> 
  summarize(lc = quantile(stab, 0.025),
            uc = quantile(stab, 0.975))

df_hsv = df |> 
  group_by(site) |> 
  summarize(var = var(hv_size),
            mean = mean(hv_size),
            sd = sd(hv_size),
            cv = sd/mean,
            stab = 1/cv) |> 
  left_join(df_mm, by = 'site') 

ggplot(df_hsv, aes(site, stab, color = site))+
  geom_point(size = 5)+
  geom_errorbar(aes(ymin = lc, ymax = uc), linewidth = 2, width = 0)+
  labs(x = 'site', y = 'Stability volume')+
  scale_color_viridis_d(option = 'turbo')+
  #scale_y_log10()+
  theme_bw()+
  theme(axis.title = element_text(size = 14), 
        axis.text = element_text(size = 14, colour = "gray0"), 
        plot.title = element_text(size = 14, hjust=0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = 'none',
        legend.title = element_text(size = 14),
        strip.text.x = element_text(size = 14),
        legend.text = element_text(size = 12))

ggsave('figs/stabhvSize.png', 
       units="in", width=8, height=6, dpi=600)

# ggsave('figs/varhvSizelog.tiff', 
#        units="in", width=8, height=6, dpi=600,compression = 'lzw')



# step 2 variable importance-----
# centroid distance----
df_ov = read_csv("sitesthrutime.csv") |> 
  select(site,y1,y2,ychange,dist_cent) |> 
  mutate(across(TT:TDR, \(x) x^2)) |> 
  pivot_longer(TT:TDR, names_to = 'axis', values_to = 'dist')

ax = unique(df_ov$axis)

for(i in 1:length(ax)){
  d = df_ov |> 
    filter(axis != ax[i]) |> 
    group_by(site,y1,y2,ychange,dist_cent) |> 
    summarise(cd = sqrt(sum(dist))) |> 
    mutate(imp = (dist_cent/cd) - 1,
           axis = ax[i])
  
  if(i == 1){
    df_imp = d
  }else{
    df_imp = bind_rows(df_imp, d)
  }
}

df_cdi = df_imp |> 
  filter(ychange == 1) |>
  group_by(site,axis) |>
  summarize(imp = mean(imp)) |>
  group_by(site) |> 
  mutate(s_imp = imp/max(imp))|> 
  mutate(site = factor(site, levels = 
                          c('JON', 'RKB', 'TWN', 'RAN', 'WHP',
                            'MAD', 'CAL', 'CRN', 'EAG', 'BLK')),
         axis = factor(axis, levels = c('TDR', 'TMA', 'SGR',
                                        'SF', 'HW', 'TT')))

y_label_formatter = function(x) {
  ifelse(x %% 1 == 0, formatC(x, format = "f", digits = 0), formatC(x, format = "f", digits = 2))
}

ggplot(df_cdi, aes(axis, s_imp, fill = site))+
  geom_col()+
  labs(x = 'Variable', y = 'Centroid distance variable importance')+
  coord_flip()+
  theme_bw()+
  facet_wrap(~site,  nrow = 2)+
  scale_y_continuous(
    breaks = c(0.0, 0.25, 0.5, 0.75, 1.0),
    limits = c(0, 1),
    labels = y_label_formatter) +
  #scale_x_discrete(labels = c(''))+
  scale_fill_viridis_d(option = 'turbo')+
  theme(axis.title = element_text(size = 14), 
        axis.text = element_text(size = 14, colour = "gray0"), 
        plot.title = element_text(size = 14, hjust=0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = 'none',
        legend.title = element_text(size = 14),
        strip.text.x = element_text(size = 14),
        legend.text = element_text(size = 12))

ggsave('figs/s_impCentDist.png', 
       units="in", width=12, height=6, dpi=600)


df_imp = readRDS('data/hvAll_vi.rds') |> 
  select(site, YEAR, hv_size, vi) |> 
  unnest_longer(vi, values_to = 'imp',indices_to = 'axis') 


df_hsi = df_imp |> 
  group_by(site,axis) |> 
  summarize(imp = mean(imp)) |> 
  group_by(site) |> 
  mutate(s_imp = imp/max(imp),
         axis = if_else(axis == 'sg_rich', 'SGR', axis))|> 
  mutate(site = factor(site, levels = 
                          c('JON', 'RKB', 'TWN', 'RAN', 'WHP',
                            'MAD', 'CAL', 'CRN', 'EAG', 'BLK')),
         axis = factor(axis, levels = c('TDR', 'TMA', 'SGR',
                                        'SF', 'HW', 'TT')))

y_label_formatter = function(x) {
  ifelse(x %% 1 == 0, formatC(x, format = "f", digits = 0), formatC(x, format = "f", digits = 2))
}

ggplot(df_hsi, aes(axis, s_imp, fill = site))+
  geom_col()+
  labs(x = 'Variable', y = 'Volume variable importance')+
  coord_flip()+
  theme_bw()+
  facet_wrap(~site,  nrow = 2)+
  scale_y_continuous(
    breaks = c(0.0, 0.25, 0.5, 0.75, 1.0),
    limits = c(0, 1),
    labels = y_label_formatter) +
  scale_fill_viridis_d(option = 'turbo')+
  theme(axis.title = element_text(size = 14), 
        axis.text = element_text(size = 14, colour = "gray0"), 
        plot.title = element_text(size = 14, hjust=0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = 'none',
        legend.title = element_text(size = 14),
        strip.text.x = element_text(size = 14),
        legend.text = element_text(size = 12))

ggsave('figs/s_imphvSize.png', 
       units="in", width=12, height=6, dpi=600)

# step 3 ----
# centroid distance
df = read_csv('hvALL.csv') %>% 
  group_by(site, YEAR) %>%
  reframe(
    var = var(hv_size),
    mean = mean(hv_size),
    sd = sd(hv_size),
    cv = sd / mean,
    stab = 1 / cv
  )

df_ov = read_csv("sitesthrutime.csv") |> 
  group_by(site, y1, y2) %>%
  reframe(
    var = var(dist_cent),
    mean = mean(dist_cent),
    sd = sd(dist_cent),
    cv = sd / mean,
    stab = 1 / cv,
    ychange = unique(ychange)) |> 
  group_by(site) |>
  nest() |> 
  mutate(m_int = map(data, \(df)lm(mean~1, data = df)),
         m_lin = map(data, \(df)lm(mean~ychange, data = df)),
         m_quad = map(data, \(df)lm(mean~ychange + I(ychange^2), data = df)),
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

ggplot(d, aes(ychange, mean, color = site))+
  geom_hline(aes(yintercept = 1), linetype = 'dashed')+
  geom_point(size = 2.5)+
  geom_smooth(data = d |> filter(model == 'Intercept'),
              method = 'lm', formula = y~1, 
              linewidth = 1, color = 'black')+
  geom_smooth(data = d |> filter(model == 'Linear'),
              method = 'lm', formula = y~x, 
              linewidth = 1, color = 'black')+
  geom_smooth(data = d |> filter(model == 'Quadratic'),
              method = 'lm', formula = y~x+I(x^2), 
              linewidth = 1, color = 'black')+
  facet_wrap(~site,  nrow = 2)+
  labs(x = 'Years between comparison', y = 'Centroid distance')+
  scale_color_viridis_d(option = 'turbo')+
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

ggsave('figs/centResp.png', 
       units="in", width=10, height=6, dpi=600)


# size ratio----
df = read_csv('hvALL.csv')

df_ov = read_csv("sitesthrutime.csv") |> 
  mutate(lsr = log(size_rat)) |> 
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

ggplot(d, aes(ychange, lsr, color = site))+
  geom_hline(aes(yintercept = 0), linetype = 'dashed')+
  geom_point(size = 2.5)+
  geom_smooth(data = d |> filter(model == 'Intercept'),
              method = 'lm', formula = y~1, 
              linewidth = 1, color = 'black')+
  geom_smooth(data = d |> filter(model == 'Linear'),
              method = 'lm', formula = y~x, 
              linewidth = 1, color = 'black')+
  geom_smooth(data = d |> filter(model == 'Quadratic'),
              method = 'lm', formula = y~x+I(x^2), 
              linewidth = 1, color = 'black')+
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


# size relative----
df = read_csv('hvAll.csv')%>% 
  group_by(site, YEAR) %>%
  reframe(
    var = var(hv_size),
    mean = mean(hv_size),
    sd = sd(hv_size),
    cv = sd / mean,
    stab = 1 / cv
  )

df_ov = read_csv("sitesthrutime.csv") |>  
  group_by(site, y1, y2) %>%
  reframe(
    var_hv1 = var(hv1_size),
    mean_hv1 = mean(hv1_size),
    sd_hv1 = sd(hv1_size),
    cv_hv1 = sd_hv1/ mean_hv1,
    var_hv2 = var(hv2_size),
    mean_hv2 = mean(hv2_size),
    sd_hv2 = sd(hv2_size),
    cv_hv2 = sd_hv2/ mean_hv2,
    stab_hv2 = 1 / cv_hv2,
    ychange = unique(ychange)) |> 
  mutate(size_ch = (mean_hv1-mean_hv2)/mean_hv1) |> 
  group_by(site) |>
  nest() |> 
  mutate(m_int = map(data, \(df)lm(size_ch~1, data = df)),
         m_lin = map(data, \(df)lm(size_ch~ychange, data = df)),
         m_quad = map(data, \(df)lm(size_ch~ychange + I(ychange^2), data = df)),
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

ggplot(d, aes(ychange, size_ch, color = site))+
  geom_hline(aes(yintercept = 1), linetype = 'dashed')+
  geom_point(size = 2.5)+
  geom_smooth(data = d |> filter(model == 'Intercept'),
              method = 'lm', formula = y~1, 
              linewidth = 1, color = 'black')+
  geom_smooth(data = d |> filter(model == 'Linear'),
              method = 'lm', formula = y~x, 
              linewidth = 1, color = 'black')+
  geom_smooth(data = d |> filter(model == 'Quadratic'),
              method = 'lm', formula = y~x+I(x^2), 
              linewidth = 1, color = 'black')+
  facet_wrap(~site,  nrow = 2)+
  labs(x = 'Years between comparison', y = 'Y2/Y1 size ratio')+theme_bw()+
  theme(axis.title = element_text(size = 14), 
        axis.text = element_text(size = 14, colour = "gray0"), 
        plot.title = element_text(size = 14, hjust=0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = 'right',
        legend.title = element_text(size = 14),
        strip.text.x = element_text(size = 14),
        legend.text = element_text(size = 12))

ggsave('figs/sizeRespLog.tiff', 
       units="in", width=13, height=9, dpi=600,compression = 'lzw')
