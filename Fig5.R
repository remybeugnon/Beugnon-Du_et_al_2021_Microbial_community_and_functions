# platform       x86_64-w64-mingw32          
# arch           x86_64                      
# os             mingw32                     
# system         x86_64, mingw32             
# language       R                           
# version.string R version 4.0.3 (2020-10-10)
# Packages:
#   - arm          1.11.2
#   - tidyverse    1.3.0
#   - lavaan       0.6-7
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

#### Packages ####
libs <- c(
  'arm','tidyverse',
  'lavaan'
)
invisible(lapply(libs, library, character.only = T))

#### > 0.2 Loading data ####
df <- read.csv(file = 'df.csv')

#### > 1. SEM fit ####
VAR = c('total.biomass', 'active.biomass', 
        'bf', 'shannon.diversity.bacteria', 'shannon.diversity.fungi', 
        'Cata', 'FG.evenness',
        'SU.range', 'SU.eff', 
        'Basal.respiration', 
        'TOC', 'CP', 'CN', 'pH', 'H2O', 'tree.species.richness')

df = apply(df[,VAR],2,rescale)

model = '
Basal.respiration ~ total.biomass + active.biomass + 
  bf + shannon.diversity.bacteria + shannon.diversity.fungi +  
  Cata  +  FG.evenness +
  SU.eff + SU.range + 
  TOC + CP + CN + pH + H2O +
  tree.species.richness

total.biomass ~~ active.biomass
total.biomass ~~ bf
total.biomass ~~ shannon.diversity.bacteria
total.biomass ~~ shannon.diversity.fungi
total.biomass ~~ Cata
total.biomass ~~ FG.evenness

active.biomass ~~ bf
active.biomass ~~ shannon.diversity.bacteria
active.biomass ~~ shannon.diversity.fungi
active.biomass ~~ Cata
active.biomass ~~ FG.evenness

bf ~~ shannon.diversity.bacteria
bf ~~ shannon.diversity.fungi
bf ~~ Cata
bf ~~ FG.evenness

shannon.diversity.bacteria ~~ shannon.diversity.fungi
shannon.diversity.bacteria ~~ Cata
shannon.diversity.bacteria ~~ FG.evenness
shannon.diversity.fungi ~~ Cata
shannon.diversity.fungi ~~ FG.evenness

Cata ~~ FG.evenness

SU.eff ~ total.biomass + active.biomass + 
         bf + shannon.diversity.bacteria + shannon.diversity.fungi +  
         Cata + FG.evenness +
         TOC + CP + CN + pH + H2O + 
         tree.species.richness
         
SU.range ~ total.biomass + active.biomass +
           bf + shannon.diversity.bacteria + shannon.diversity.fungi +  
           Cata  + FG.evenness +
           TOC + CP + CN + pH + H2O + 
           tree.species.richness

SU.range ~~ SU.eff

Cata ~ TOC + CP + CN + pH + H2O + tree.species.richness
FG.evenness ~ TOC + CP + CN + pH + H2O + tree.species.richness

shannon.diversity.fungi ~ TOC + CP + CN + pH + H2O + tree.species.richness
shannon.diversity.bacteria ~ TOC + CP + CN + pH + H2O + tree.species.richness
bf ~ TOC + CP + CN + pH + H2O + tree.species.richness

active.biomass ~ TOC + CP + CN + pH + H2O + tree.species.richness
total.biomass ~ TOC + CP + CN + pH + H2O + tree.species.richness

tree.species.richness ~~ CP
tree.species.richness ~~ CN
tree.species.richness ~~ pH
tree.species.richness ~~ H2O
tree.species.richness ~~ TOC

TOC ~~ CP
TOC ~~ CN
TOC ~~ pH
TOC ~~ H2O

CP ~~ CN
CP ~~ pH
CP ~~ H2O

CN ~~ pH
CN ~~ H2O

pH ~~ H2O
'

fit = sem(model,df, fixed.x = F)

#### > 2. Plot R.squared ####
r = c("Microbial resp.", 'SU eff.', 'SU range', 'Cata','FG eve.' , 'Bac. div.', 'Fung. div.', 'B:F', 'Active biomass', 'Biomass')
R2 = inspect(fit, 'r2') %>% round(digits = 2) %>% data.frame()
R2$labels = r %>% factor(., levels = r)
R2$. = (R2$. * 100)

p.r.squaed =
  ggplot(data = R2, 
         aes(y = ., 
             x = labels, 
             fill = labels)) +
  geom_hline(aes(yintercept = 0), color = 'gray', linetype = 2) +
  geom_hline(aes(yintercept = 25), color = 'gray', linetype = 2) +
  geom_hline(aes(yintercept = 50), color = 'gray', linetype = 2) +
  geom_bar(stat = 'identity', width = 0.5) +
  labs(x = '', 
       y = '') +
  coord_flip() +
  scale_fill_manual(breaks = R2$labels,
                     values = c("Microbial resp." = '#52A55D', 
                                'SU eff.' = '#D4C86A', 'SU range' = '#D4C86A', 
                                'Cata' = '#B45A81', 'FG eve.' = '#B45A81',
                                'Bac. div.' = '#675091', 'Fung. div.' = '#675091', 'B:F' = '#675091', 
                                'Active biomass' = '#675091', 'Biomass' = '#675091')) +
  scale_y_continuous(breaks = c(0,25,50,75))+
  theme(panel.grid = element_blank(), 
          plot.background = element_rect(fill = 'white', color = 'white'), 
          panel.background = element_rect(fill = 'gray95', color = 'white'), 
          legend.position = 'none', 
          axis.text.x = element_text(size = 9),
          axis.text.y = element_text(size = 9), plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"))

ggsave(p.r.squaed, filename = 'Fig5-r2.png', height = 11, width = 4, units = 'cm')
#### Figure 5 ####
# The final figure have been made manually using Inkscape