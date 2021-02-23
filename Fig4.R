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

#### > 0.1 Packages ####
libs <- c(
  'arm','readxl',
  'tidyverse',
  'lavaan'
)
invisible(lapply(libs, library, character.only = T))

#### > 0.2 Loading data ####
df <- read.csv(file = 'Data/df-models-vf2.csv')

#### > 1. Statistical analyses  ####
Exp = c('total.biomass', 'active.biomass', 'bf', 'shannon.diversity.bacteria', 'shannon.diversity.fungi')
Resp = c('microbial.respiration' , 'SU.range', 'SU.eff', 'Cata', 'FG.evenness')

df = apply(df[,c(Exp, Resp)],2,rescale)

model = '
microbial.respiration ~ total.biomass + active.biomass + 
                        bf + shannon.diversity.bacteria + shannon.diversity.fungi +  
                        Cata + FG.evenness +  SU.eff + SU.range 

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
         Cata + FG.evenness
SU.range ~ total.biomass + active.biomass + 
           bf + shannon.diversity.bacteria + shannon.diversity.fungi +  
           Cata + FG.evenness
           
SU.range ~~ SU.eff
'

fit = sem(model,df)
params <- lavaan::standardizedSolution(fit)

# Colors
non.sign = 'gray80'

# Node properties
r2 = inspect(fit, 'r2')

# Add R.suared values to node names
shape.plot = read_xlsx('Plot_SEM_Shape.xlsx', sheet = 'Nodes')
shape.plot$label[shape.plot$node == "microbial.respiration"] = paste0('Microbial respiration (',(r2['microbial.respiration']*100)%>%round(1),'%)' )
shape.plot$label[shape.plot$node == "SU.range"] = paste0('SU range\n(',(r2['SU.range']*100)%>%round(1),'%)' )
shape.plot$label[shape.plot$node == "SU.eff"] = paste0('SU efficiency\n(',(r2['SU.eff']*100)%>%round(1),'%)' )
shape.plot$label[shape.plot$node == "Cata"] = 'Cata'
shape.plot$label[shape.plot$node == "FG.evenness"] = 'FG eve.'
shape.plot$label[shape.plot$node == "shannon.diversity.bacteria"] = 'Bac.\ndiv.'
shape.plot$label[shape.plot$node == "shannon.diversity.fungi"] = 'Fung.\ndiv.'
shape.plot$label[shape.plot$node == "bf"] = 'B:F'
shape.plot$label[shape.plot$node == "active.biomass"] = 'Active\nbiomass'
shape.plot$label[shape.plot$node == "total.biomass"] = 'Total\nbiomass'

# Extracting edges from fit
param_edges <- params %>%
  filter(op %in% c("=~", "~", "~~"), lhs != rhs) %>%
  transmute(to = lhs,
            from = rhs,
            val = est.std,
            pval = pvalue,
            type = dplyr::case_when(
              op == "=~" ~ "loading",
              op == "~"  ~ "regression",
              op == "~~" ~ "correlation",
              TRUE ~ NA_character_))

# Addition of positions
edge.plot = read_xlsx('Plot_SEM_Shape.xlsx', sheet = 'Edge') 
edge.plot[,3:6] = round(edge.plot[,3:6],1)
edge.plot = left_join(edge.plot, param_edges, by = c('from','to'))

# Significance
edge.plot$si = '' 
edge.plot$si[edge.plot$pval <= 0.05 & edge.plot$pval > 0.01] = ' *'
edge.plot$si[edge.plot$pval <= 0.01 & edge.plot$pval > 0.001] = ' **'
edge.plot$si[edge.plot$pval <= 0.001] = ' ***'

# calculation of standard positions
edge.plot$label.x = edge.plot$start.x - 0.5 * (edge.plot$start.x - edge.plot$stop.x) 
edge.plot$label.y = edge.plot$start.y - 0.5 * (edge.plot$start.y - edge.plot$stop.y)
edge.plot$angle.text = 0

# Manual corrections
edge.plot$angle.text[edge.plot$from == 'total.biomass' & edge.plot$to == 'SU.eff'] = 90
edge.plot$label.x[edge.plot$from == 'total.biomass' & edge.plot$to == 'SU.eff'] = 10.9
edge.plot$label.y[edge.plot$from == 'total.biomass' & edge.plot$to == 'SU.eff'] = 10.9
edge.plot$angle.text[edge.plot$from == 'FG.evenness' & edge.plot$to == 'bf'] = 64.5
edge.plot$label.y[edge.plot$from == 'FG.evenness' & edge.plot$to == 'bf'] = 15.6  
edge.plot$label.x[edge.plot$from == 'FG.evenness' & edge.plot$to == 'bf'] = 15

edge.plot$angle.text[edge.plot$from == 'active.biomass' & edge.plot$to == 'microbial.respiration'] = 90
edge.plot$label.x[edge.plot$from == 'active.biomass' & edge.plot$to == 'microbial.respiration'] = 8.2
edge.plot$label.y[edge.plot$from == 'active.biomass' & edge.plot$to == 'microbial.respiration'] = 10.9
edge.plot$angle.text[edge.plot$from == 'active.biomass' & edge.plot$to == 'SU.eff'] = 90
edge.plot$label.x[edge.plot$from == 'active.biomass' & edge.plot$to == 'SU.eff'] = 9.7
edge.plot$label.y[edge.plot$from == 'active.biomass' & edge.plot$to == 'SU.eff'] = 10.9

edge.plot$label.y[edge.plot$from == 'shannon.diversity.fungi' & edge.plot$to == 'active.biomass'] = 13.45
edge.plot$label.x[edge.plot$from == 'shannon.diversity.fungi' & edge.plot$to == 'active.biomass'] = 15.8
edge.plot$angle.text[edge.plot$from == 'shannon.diversity.fungi' & edge.plot$to == 'total.biomass'] = -10
edge.plot$label.x[edge.plot$from == 'shannon.diversity.fungi' & edge.plot$to == 'total.biomass'] = 17.8
edge.plot$label.y[edge.plot$from == 'shannon.diversity.fungi' & edge.plot$to == 'total.biomass'] = 16.6
edge.plot$angle.text[edge.plot$from == 'shannon.diversity.fungi' & edge.plot$to == 'microbial.respiration'] = -90
edge.plot$label.x[edge.plot$from == 'shannon.diversity.fungi' & edge.plot$to == 'microbial.respiration'] = 21.8
edge.plot$label.y[edge.plot$from == 'shannon.diversity.fungi' & edge.plot$to == 'microbial.respiration'] = 6

edge.plot$angle.text[edge.plot$from == 'bf' & edge.plot$to == 'microbial.respiration'] = -90
edge.plot$label.x[edge.plot$from == 'bf' & edge.plot$to == 'microbial.respiration'] = 22.8
edge.plot$label.y[edge.plot$to == 'bf' & edge.plot$from == 'Cata'] = 14.5
edge.plot$label.x[edge.plot$to == 'bf' & edge.plot$from == 'Cata'] = 18.3
edge.plot$angle.text[edge.plot$to == 'bf' & edge.plot$from == 'Cata'] = -90
edge.plot$angle.text[edge.plot$from == 'bf' & edge.plot$to == 'active.biomass'] = 13.8
edge.plot$label.x[edge.plot$from == 'bf' & edge.plot$to == 'active.biomass'] = 14.3
edge.plot$label.y[edge.plot$from == 'bf' & edge.plot$to == 'active.biomass'] = 16.35
edge.plot$label.y[edge.plot$from == 'bf' & edge.plot$to == 'total.biomass'] = 18.3

edge.plot$angle.text[edge.plot$from == 'FG.evenness' & edge.plot$to == 'total.biomass'] = 90
edge.plot$label.x[edge.plot$from == 'FG.evenness' & edge.plot$to == 'total.biomass'] = 12.2
edge.plot$label.y[edge.plot$from == 'FG.evenness' & edge.plot$to == 'total.biomass'] = 14.9
edge.plot$angle.text[edge.plot$from == 'FG.evenness' & edge.plot$to == 'SU.eff'] = 90
edge.plot$label.x[edge.plot$from == 'FG.evenness' & edge.plot$to == 'SU.eff'] = 12.2
edge.plot$angle.text[edge.plot$from == 'FG.evenness' & edge.plot$to == 'microbial.respiration'] = 90
edge.plot$label.x[edge.plot$from == 'FG.evenness' & edge.plot$to == 'microbial.respiration'] = 13
edge.plot$label.y[edge.plot$from == 'FG.evenness' & edge.plot$to == 'microbial.respiration'] = 4.5

edge.plot$angle.text[edge.plot$from == 'Cata' & edge.plot$to == 'microbial.respiration'] = -90
edge.plot$label.x[edge.plot$from == 'Cata' & edge.plot$to == 'microbial.respiration'] = 17
edge.plot$label.y[edge.plot$from == 'Cata' & edge.plot$to == 'microbial.respiration'] = 7.9
edge.plot$angle.text[edge.plot$from == 'Cata' & edge.plot$to == 'SU.range'] = -90
edge.plot$label.x[edge.plot$from == 'Cata' & edge.plot$to == 'SU.range'] = 18.5
edge.plot$angle.text[edge.plot$from == 'Cata' & edge.plot$to == 'SU.eff'] = 30
edge.plot$label.x[edge.plot$from == 'Cata' & edge.plot$to == 'SU.eff'] = 14.55
edge.plot$label.y[edge.plot$from == 'Cata' & edge.plot$to == 'SU.eff'] = 8.2
edge.plot$label.y[edge.plot$to == 'Cata' & edge.plot$from == 'FG.evenness'] = 10.5
edge.plot$angle.text[edge.plot$from == 'Cata' & edge.plot$to == 'total.biomass'] = -55
edge.plot$label.y[edge.plot$from == 'Cata' & edge.plot$to == 'total.biomass'] = 12.25
edge.plot$label.x[edge.plot$from == 'Cata' & edge.plot$to == 'total.biomass'] = 16.1


edge.plot$angle.text[edge.plot$from == 'SU.eff' & edge.plot$to == 'microbial.respiration'] = 90
edge.plot$label.x[edge.plot$from == 'SU.eff' & edge.plot$to == 'microbial.respiration'] = 12.2
edge.plot$angle.text[edge.plot$from == 'SU.range' & edge.plot$to == 'microbial.respiration'] = -90
edge.plot$label.x[edge.plot$from == 'SU.range' & edge.plot$to == 'microbial.respiration'] = 17.8
edge.plot$label.y[edge.plot$from == 'SU.range' & edge.plot$to == 'SU.eff'] = 6.3

# Plot
p =
  ggplot(shape.plot, aes(x = center.x, y = center.y)) +
  # Background rectangles
  geom_rect(aes(xmin = 0, xmax = 24 , ymin = 8, ymax = 21), fill = '#C0C2EE') + 
  geom_rect(aes(xmin = 0, xmax = 24 , ymin = 4.5, ymax = 7.5), fill = '#FFF0C5') + 
  geom_rect(aes(xmin = 0, xmax = 24 , ymin = .5, ymax = 3.5), fill = '#EDFCC3') +
  geom_rect(aes(xmin = 6.5, xmax = 14.5 , ymin = 13.65, ymax = 20.7), color = '#414263', alpha = 0.5, fill = NA, lty = 5) +
  geom_rect(aes(xmin = 15.5, xmax = 23.5 , ymin = 13.65, ymax = 20.7), color = '#414263', alpha = 0.5, fill = NA, lty = 5) + 
  geom_rect(aes(xmin = 11.4, xmax = 23.5 , ymin = 8.5, ymax = 13), color = '#414263', alpha = 0.5, fill = NA, lty = 5) +
  # Nodes 
  geom_rect(aes(xmin = (center.x - (width / 2)),
                xmax = (center.x + (width / 2)),
                ymin = (center.y - (height / 2)),
                ymax = (center.y + (height / 2))),
            color = 'white', fill = 'white', alpha = 0.5) +
  # Non-significant 
  ## Correlation
  geom_segment(data = edge.plot %>% filter(curve != 'y' & type == 'correlation' & pval > 0.05), 
               aes(x = start.x, y = start.y, 
                   xend = stop.x, yend = stop.y), 
               size = .4,
               arrow = arrow(length = unit(0.01, "npc"), ends = 'both', angle=30), 
               lineend = 'butt',  lty = 2,
               color = non.sign) +  
  ## Correlation curve
  geom_curve(data = edge.plot %>% filter(curve == 'y' & type == 'correlation' & pval > 0.05), 
             aes(x = start.x, y = start.y, 
                 xend = stop.x, yend = stop.y), 
             size = .4, 
             lineend = 'butt',  lty = 2,
             arrow = arrow(length = unit(0.01, "npc"), ends = 'both', angle=30), 
             color = non.sign,
             curvature = -.2) +
  ## Regression
  geom_segment(data = edge.plot %>% filter(curve != 'y' & type == 'regression' & pval > 0.05), 
               aes(x = start.x, y = start.y, 
                   xend = stop.x, yend = stop.y),  
               size = .4,
               lineend = 'butt', linejoin = 'mitre', lty = 2,
               arrow = arrow(length = unit(0.01, "npc"), angle=30), color = non.sign) +  
  # Regression curve
  geom_curve(data = edge.plot %>% filter(curve == 'y' & type == 'regression' & pval > 0.05), 
             aes(x = start.x, y = start.y, 
                 xend = stop.x, yend = stop.y),  
             size = .4,
             lineend = 'butt', lty = 2,
             arrow = arrow(length = unit(0.01, "npc"), angle=30), color = non.sign,
             curvature = -.2) +  
  # Significant
  ## Regression
  geom_segment(data = edge.plot %>% filter(curve != 'y' & type == 'regression' & pval <= 0.05), 
               aes(x = start.x, y = start.y, 
                   xend = stop.x, yend = stop.y,
                   size = abs(val), color = val), 
               lineend = 'butt', linejoin = 'mitre',
               arrow = arrow(length = unit(0.01, "npc"), angle=30)) +
  
  
  ## Regression curve
  geom_curve(data = edge.plot %>% filter(curve == 'y' & type == 'regression' & pval <= 0.05), 
             aes(x = start.x, y = start.y, 
                 xend = stop.x, yend = stop.y,
                 size = abs(val), color = val), 
             lineend = 'butt',
             arrow = arrow(length = unit(0.01, "npc"), angle=30),
             curvature = -.2) +
  # Correlation
  geom_segment(data = edge.plot %>% filter(curve != 'y' & type == 'correlation' & pval <= 0.05), 
               aes(x = start.x, y = start.y, 
                   xend = stop.x, yend = stop.y,
                   size = abs(val), color = val), 
               lineend = 'butt', 
               arrow = arrow(length = unit(0.01, "npc"), ends = 'both', angle=30)) +
  # Correlation curve
  geom_curve(data = edge.plot %>% filter(curve == 'y' & type == 'correlation' & pval <= 0.05), 
             aes(x = start.x, y = start.y, 
                 xend = stop.x, yend = stop.y,
                 size = abs(val), color = val), 
             lineend = 'butt', 
             arrow = arrow(length = unit(0.01, "npc"), ends = 'both', angle=30),
             curvature = -.2) +
  # Node text
  annotate(geom = 'text', x = shape.plot$center.x, y = shape.plot$center.y, label = shape.plot$label, size = 5.5) + 
  
  # Groups text
  annotate(geom = 'text', x = 10.5, y = 20, label = "Microbial biomass", size = 7) + 
  annotate(geom = 'text', x = 19.75, y = 20, label = "Taxonomic profile", size = 7) +
  annotate(geom = 'text', x = 21, y = 11, label = "Functional\nprofile", size = 7) + 
  annotate(geom = 'text', x = 3, y = 2, label = "Ecosystem\nfunction", size = 8) + 
  annotate(geom = 'text', x = 3, y = 6, label = "Physiological\npotential", size = 8) +
  annotate(geom = 'text', x = 3, y = 14.75, label = "Microbial\ncommunity", size = 8) +
  
  # Estimates
  geom_text(data = edge.plot %>% filter(pval<=0.05), 
            aes(x = label.x, 
                y = label.y, 
                angle = angle.text,
                label = paste0(val %>% round(3), si)),
            size = 4) +
  
  # Nice looking
  labs(x = '', y = '') + 
  scale_colour_gradient2(guide = FALSE, low = "#DD2E29", mid = "darkgray", high = "#20B120") +
  scale_size(guide = FALSE, range = c(0.5,3)) +
  lims(x = c(0,24.5), y = c(0,21.5)) +
  theme(plot.margin = margin(0,0,0,0), 
        title = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        plot.background = element_rect(color = 'white', fill = 'white'),
        panel.background =  element_rect(color = 'white', fill = 'white')
  )

ggsave(p, filename = 'Figure4.png', width = 25.5, height = 22.5, units = 'cm')