# platform       x86_64-w64-mingw32          
# arch           x86_64                      
# os             mingw32                     
# system         x86_64, mingw32             
# language       R                           
# version.string R version 4.0.3 (2020-10-10)
# Packages:
#   - arm          1.11.2
#   - tidyverse    1.3.0
#   - readxl       1.3.1
#   - ggpubr       0.4.0 
# Attached tidyverse:
#       - ggplot2 3.3.3
#       - purrr   0.3.4
#       - tibble  3.0.6 
#       - dplyr   1.0.4
#       - tidyr   1.1.2
#       - stringr 1.4.0
#       - readr   1.4.0
#       - forcats 0.5.1

rm(list = ls())

#### >> 0.1 Packages ####
libs <- c(
  'arm','tidyverse', 'readxl', 'ggpubr'
)
invisible(lapply(libs, library, character.only = T))

#### >> 0.2 Loading data ####
df <- read.csv('Data/df.csv')

#### > 1. Models ####
Exp.1 = c('total.biomass', 'active.biomass', 'bf', 'Shannon.div.bac', 'Shannon.div.fung')
Resp.1 = c('Basal.respiration' , 
           'SU.range', 'SU.eff',
           'ab.fg.C', 'r.fg', 'fg.div', 'Cata', 'r.fg.r', 'FG.evenness')

# Models
MOD = list(
  mod.bio = lm(data = df, formula = total.biomass ~ tree.species.richness),
  mod.act.bio = lm(data = df, formula = active.biomass ~ tree.species.richness),
  mod.bf = lm(data = df, formula = bf ~ tree.species.richness),
  mod.bac.div = lm(data = df, formula = Shannon.div.bac ~ tree.species.richness),
  mod.fung.div = lm(data = df, formula = Shannon.div.fung ~ tree.species.richness),
  mod.cata = lm(data = df, formula = Cata ~ tree.species.richness),
  mod.fg.eve = lm(data = df, formula = FG.evenness ~ tree.species.richness),
  mod.SUeff = lm(data = df, formula = SU.eff ~ tree.species.richness),
  mod.SUrange = lm(data = df, formula = SU.range ~ tree.species.richness),
  mod.basal = lm(data = df, formula = Basal.respiration ~ tree.species.richness)
)

#### > 2. Plot results ####
df.coef = bind_rows(
  data.frame(resp = '    ', summary(MOD[['mod.basal']])$coefficients),  
  data.frame(resp = 'SU range', summary(MOD[['mod.SUrange']])$coefficients),
  data.frame(resp = 'SU\nefficiency', summary(MOD[['mod.SUeff']])$coefficients),
  data.frame(resp = 'FG div.', summary(MOD[['mod.fg.eve']])$coefficients),
  data.frame(resp = 'Cata.', summary(MOD[['mod.cata']])$coefficients),
  data.frame(resp = 'Fung. div.', summary(MOD[['mod.fung.div']])$coefficients),
  data.frame(resp = 'Bac. div.', summary(MOD[['mod.bac.div']])$coefficients),
  data.frame(resp = 'B:F', summary(MOD[['mod.bf']])$coefficients),
  data.frame(resp = 'Active\nbiomass', summary(MOD[['mod.act.bio']])$coefficients),
  data.frame(resp = 'Total\nbiomass', summary(MOD[['mod.bio']])$coefficients)
  )

df.coef = df.coef[grep("tree.species.richness", rownames(df.coef)),]
df.coef$resp = df.coef$resp %>% factor(., levels = df.coef$resp)

df.signi = df.coef %>% 
  filter(Pr...t..<0.1)
df.signi$sn = '    **'
df.signi$sn[df.signi$Pr...t..<0.001] = '    ***'
df.signi$sn[df.signi$Pr...t..>0.01] = '    *'
df.signi$sn[df.signi$Pr...t..>0.05] = '    .'

#### >> 2.1 Plot estimates ####
p.est = 
  ggplot(data = df.coef, aes(x = resp , y = Estimate)) +
  geom_point(data = df.coef, aes(x = resp , y = Estimate)) + 
  geom_errorbar(data = df.coef, aes(ymin = (Estimate - 1.96 * Std..Error), ymax = (Estimate + 1.96 * Std..Error)), width = .05) +
  geom_hline(yintercept = 0) +
  geom_text(data = df.signi, aes(x = resp, y = (Estimate + 1.96 * Std..Error), label = sn), size = 5) +
  ylim(-.05,.05) +
  coord_flip() +
  labs(x = '', y = 'Tree species richness\neffect' ) +
  theme_classic()


#### >> 2.2 Plot relationships ####
#### >>> 2.2.1 Plot microbial respiration ####
gd <- df %>% 
  group_by(tree.species.richness) %>% 
  summarise(
    mean = mean(Basal.respiration),
    int.p = mean(Basal.respiration) + 1.96 * sd(Basal.respiration) / sqrt(n()),
    int.n = mean(Basal.respiration) - 1.96 * sd(Basal.respiration) / sqrt(n())  )

p.basal =
  ggplot(data = df, aes(x = tree.species.richness, y = Basal.respiration)) +
    geom_jitter(data = df, aes(x = tree.species.richness, y = Basal.respiration), size = 1) +
    geom_smooth(data = df, aes(x = tree.species.richness, y = Basal.respiration), method = 'lm', color = 'gray') +
    geom_point(data = gd, aes(x = tree.species.richness, y = mean), size = 3) + 
    geom_errorbar(data = gd, aes(x = tree.species.richness, y = mean, ymin = int.n, ymax = int.p)) +
    scale_x_continuous(breaks = c(1,4,8,16,24)) + 
    labs(x = 'p.value = 0.064.', y = bquote(atop("Microbial respiration","["~mu~mol[O2]~' '~h^{-1}~' '~g^{-1}~']'))) + 
    theme_classic()

#### >>> 2.2.2 Plot SU efficiency ####
gd <- df %>% 
  group_by(tree.species.richness) %>% 
  summarise(
    mean = mean(SU.eff),
    int.p = mean(SU.eff) + 1.96 * sd(SU.eff) / sqrt(n()),
    int.n = mean(SU.eff) - 1.96 * sd(SU.eff) / sqrt(n()))

p.SUeff =
  ggplot(data = df, aes(x = tree.species.richness, y = SU.eff)) +
  geom_jitter(data = df, aes(x = tree.species.richness, y = SU.eff), size = 1) +
  geom_smooth(data = df, aes(x = tree.species.richness, y = SU.eff), method = 'lm', color = 'gray') +
  geom_point(data = gd, aes(x = tree.species.richness, y = mean), size = 3) + 
  geom_errorbar(data = gd, aes(x = tree.species.richness, y = mean, ymin = int.n, ymax = int.p)) +
  scale_x_continuous(breaks = c(1,4,8,16,24)) + 
  labs(x = 'p.value =  0.001**', y = "SU efficiency") +
  theme_classic()

#### >>> 2.2.3 Plot Bacterial diversity ####
gd <- df %>% 
  group_by(tree.species.richness) %>% 
  summarise(
    mean = mean(Shannon.div.bac),
    int.p = mean(Shannon.div.bac) + 1.96 * sd(Shannon.div.bac) / sqrt(n()),
    int.n = mean(Shannon.div.bac) - 1.96 * sd(Shannon.div.bac) / sqrt(n()))

p.bac.div =
  ggplot(data = df, aes(x = tree.species.richness, y = Shannon.div.bac)) +
  geom_jitter(data = df, aes(x = tree.species.richness, y = Shannon.div.bac), size = 1) +
  geom_smooth(data = df, aes(x = tree.species.richness, y = Shannon.div.bac), method = 'lm', color = 'gray') +
  geom_point(data = gd, aes(x = tree.species.richness, y = mean), size = 3) + 
  geom_errorbar(data = gd, aes(x = tree.species.richness, y = mean, ymin = int.n, ymax = int.p)) +
  scale_x_continuous(breaks = c(1,4,8,16,24)) + 
  labs(x = 'p.value =  0.011*', y = "Bacteria diversity") +
  theme_classic()

#### >>> 2.2.4 Plot total biomass ####
gd <- df %>% 
  group_by(tree.species.richness) %>% 
  summarise(
    mean = mean(total.biomass),
    int.p = mean(total.biomass) + 1.96 * sd(total.biomass) / sqrt(n()),
    int.n = mean(total.biomass) - 1.96 * sd(total.biomass) / sqrt(n()))

p.bio =
  ggplot(data = df, aes(x = tree.species.richness, y = total.biomass/1000)) +
  geom_jitter(data = df, aes(x = tree.species.richness, y = total.biomass/1000), size = 1) +
  geom_smooth(data = df, aes(x = tree.species.richness, y = total.biomass/1000), method = 'lm', color = 'gray') +
  geom_point(data = gd, aes(x = tree.species.richness, y = mean/1000), size = 3) + 
  geom_errorbar(data = gd, aes(x = tree.species.richness, y = mean/1000, ymin = int.n/1000, ymax = int.p/1000)) +
  scale_x_continuous(breaks = c(1,4,8,16,24)) + 
  labs(x = 'p.value =  0.003**', y = bquote(atop("Total biomass","["~10^3~'.'~mu~g[C]~' '~g^{-1}~']'))) +
  theme_classic()
  
#### >> 2.3 Final plot ####
p2 = ggarrange(p.basal, p.SUeff, p.bac.div, p.bio, align = 'hv') %>%
  annotate_figure(bottom = "   ")

p = ggarrange(p.est, p2, widths = c(0.4, 0.6))
p

ggsave(p, 
       filename = 'Fig2-Diversity-effect.tiff',
       width = 20,
       height = 10, 
       units = 'cm'
       )

#### Figure 2 ####
# The final figure combining all the plot have been made manually using Inkscape