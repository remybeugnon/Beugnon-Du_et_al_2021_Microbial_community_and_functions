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
#   - eulerr       6.1.0
#   - ggplotify    0.0.5 
#   - cowplot      1.1.1
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
  'arm','tidyverse', 'readxl',
  'eulerr', 'ggplotify', 'cowplot'
)
invisible(lapply(libs, library, character.only = T))

#### >> 0.2 Loading data ####
df <- read.csv('df.csv')

#### > 1 Correlation between microbial facets ####
d.cor = df %>% 
  select("Total\nbiomass" = 'total.biomass', 'Active\nbiomass'='active.biomass',
         'B:F' = 'bf' , 'Bacteria\ndiversity' = 'shannon.diversity.bacteria', 'Fungi\ndiversity' = 'shannon.diversity.fungi',
         'Cata' = 'Cata', 'FG\nevenness'='FG.evenness')

m.cor = cor(d.cor)
test.cor = cor.mtest(d.cor)

col1 <- colorRampPalette(c("Darkred", "white", "Darkblue"))

# Plot
tiff('Fig3_cortable.tiff', height = 20, width = 20, units = 'cm', res = 600)

corrplot::corrplot(m.cor, method = 'color',
                   p.mat = test.cor$p, insig = "label_sig",type = "lower",
                   sig.level = c(.001, .01, .05), 
                   pch.cex = .9, pch.col = "black",
                   diag = F, 
                   tl.col = 'black', tl.srt = 0, tl.offset = 2,
                   col = col1(1000))

dev.off()

#### > 2 Microbial facets effects on microbial functions ####
#### >> 2.1 Model selection ####
#### >>> 2.1.1 Models ####
Exp = c('total.biomass', 'active.biomass', 'bf',
        'shannon.diversity.bacteria', 'shannon.diversity.fungi',
        'cata','fg.evenness')

Resp = c('microbial.respiration' , 
         'SU.range', 'SU.eff'
         )

# Scaling variables
df[,c(Exp, Resp)] = apply(df[,c(Exp, Resp)],2,rescale)

# Models 
MOD.RESP = list()
for(i in 1:length(Resp)){
  mod = paste(Resp[i], ' ~ ', paste(Exp, collapse =  ' + ' ) ) %>% 
    formula %>%
    lm(data = df, formula = .) %>%
    step(direction = 'both', trace = F)
  MOD.RESP[[i]] = list(Resp[i], mod)
}

#### >>> 2.1.2 Extraction estimates ####
exp = rep(Exp, length(Resp))
df.plot = data.frame(
  exp = exp[order(exp)],
  resp = rep(Resp, length(unique(Exp))), 
  estimate = NA,
  pv = NA,
  r2 = NA,
  stringsAsFactors = F
)

for(i in 1:length(Resp)){
  resp = MOD.RESP[[i]][[1]]
  mod =  MOD.RESP[[i]][[2]] %>% summary
  df.plot$r2[df.plot$resp == resp] = mod$adj.r.squared * 100
  mod.coeff = mod$coefficients %>% data.frame()
  exp = mod.coeff %>% row.names() %>% .[.!= '(Intercept)']
  for(j in exp){
    df.plot$estimate[df.plot$exp == j & df.plot$resp == resp] = mod.coeff$Estimate[rownames(mod.coeff)==j] %>% as.numeric()
    df.plot$pv[df.plot$exp == j & df.plot$resp == resp] = mod.coeff$Pr...t..[rownames(mod.coeff)==j] %>% as.numeric()
  }
}

#### >> 2.2. Variance partitioning ####
#### >>> 2.2.1 Variance partitioning model ####
part.resp = varpart(Y = df$microbial.respiration, 
                    df%>%dplyr::select(active.biomass), 
                    df%>%select(bf,shannon.diversity.bacteria,shannon.diversity.fungi), 
                    df%>%select(Cata)
)

part.SUeff = varpart(Y = df$shannon.microresp, 
                     df%>%select(ab.mic.com,active.biomass), 
                     df%>%select(shannon.diversity.bacteria,shannon.diversity.fungi), 
                     df%>%select(FG.evenness)
)

part.SUrange = varpart(Y = df$ala.oxa, 
                       df%>%select(active.biomass), 
                       df%>%select(FG.evenness)
)

#### >>> 2.2.2 Variance partitioning plots ####
# Microbial respiration
a.resp = c(part.resp$part$fract[1:3 ,'Adj.R.square'], part.resp$part$indfract[4:7 ,'Adj.R.square']) %>% unlist
a.resp[a.resp<0] = 0
names(a.resp) = NULL
names(a.resp) = c("A", "B", "C", "A&B", "B&C", "A&C", "A&B&C")
VennDiag <- euler(a.resp)
p.resp <-  plot(VennDiag, counts = TRUE, 
                font = 1, cex = .5, alpha = 0.5,
                fill=c("#403075", "#923158", "#256F5B"), 
                color = 'black',
                labels = c('B', 'T', 'F'))

# Plot
tiff('Fig3_varpart-resp.tiff', height = 5, width = 5, units = 'cm', res = 600)
p.resp 
dev.off()

# Substrate use efficiency
a.SUeff = c(part.SUeff$part$fract[1:3 ,'Adj.R.square'], part.SUeff$part$indfract[4:7 ,'Adj.R.square']) %>% unlist
a.SUeff[a.SUeff<0] = 0
names(a.SUeff) = NULL
names(a.SUeff) = c("A", "B", "C", "A&B", "B&C", "A&C", "A&B&C")
VennDiag <- euler(a.SUeff)
p.SUeff <-  plot(VennDiag, counts = TRUE, 
                 font = 1, cex = .5, alpha = 0.5,
                 fill=c("#403075", "#923158", "#256F5B"), 
                 color = 'black',
                 labels = c('B', 'T', 'F'))

tiff('Fig3_varpart-SUeff.tiff', height = 5, width = 5, units = 'cm', res = 600)
p.SUeff 
dev.off()

# Substrate use range
a.SUrange = c(part.SUrange$part$fract[1:2 ,'Adj.R.squared'], 0) %>% unlist
a.SUrange[a.SUrange<0] = 0
names(a.SUrange) = NULL
names(a.SUrange) = c("A", "B", "A&B")
VennDiag <- euler(a.SUrange)
p.SUrange <-  plot(VennDiag, counts = TRUE, 
                   font = 1, cex = 2, alpha = 0.5,
                   fill=c("#403075", "#256F5B"), 
                   color = 'black',
                   labels = c('B', 'F'))

tiff('Fig3_varpart-SUrange.tiff', height = 5, width = 5, units = 'cm', res = 600)
p.SUrange 
dev.off()

#### >>> 2.2.3 Plot model R2 ####
p.resp = ggplot(data = NULL) +
  geom_bar(aes(y = 1, x = 1), stat = 'identity', color = 'black', fill = 'white') +
  geom_bar(aes(y = part.resp$part$fract[7,'Adj.R.square'], x = 1), stat = 'identity', color = 'black', fill = 'black') +
  coord_flip(ylim = c(0,1), xlim = c(.2,1.5), clip = 'off') +
  theme_classic() +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(), 
        panel.border = element_blank()) + 
  annotate('line', y = c(0,1), x = c(.45,.45)) +  
  annotate('line', y = c(0,0), x = c(.45,.4)) +
  annotate('line', y = c(.25,.25), x = c(.45,.4)) +  
  annotate('line', y = c(.5,.5), x = c(.45,.4)) +
  annotate('line', y = c(.75,.75), x = c(.45,.4)) +
  annotate('line', y = c(1,1), x = c(.45,.4)) +  
  annotate('text', y = 0.05, x = .25, label = '0%', size = 10) + 
  annotate('text', y = .5, x = .25, label = '50%', size = 10) + 
  annotate('text', y = .95, x = .25, label = '100%', size = 10) + 
  annotate('text', y = .76, x = 1, size = 10,
           label = bquote(.(round(part.resp$part$fract[7,'Adj.R.square'], digits = 2) * 100) ~ "%"))

p.SUeff = ggplot(data = NULL) +
  geom_bar(aes(y = 1, x = 1), stat = 'identity', color = 'black', fill = 'white') +
  geom_bar(aes(y = part.SUeff$part$fract[7,'Adj.R.square'], x = 1), stat = 'identity', color = 'black', fill = 'black') +
  coord_flip(ylim = c(0,1), xlim = c(.2,1.5), clip = 'off') +
  theme_classic() +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(), 
        panel.border = element_blank()) + 
  annotate('line', y = c(0,1), x = c(.45,.45)) +  
  annotate('line', y = c(0,0), x = c(.45,.4)) +
  annotate('line', y = c(.25,.25), x = c(.45,.4)) +  
  annotate('line', y = c(.5,.5), x = c(.45,.4)) +
  annotate('line', y = c(.75,.75), x = c(.45,.4)) +
  annotate('line', y = c(1,1), x = c(.45,.4)) +  
  annotate('text', y = 0.05, x = .25, label = '0%', size = 10) + 
  annotate('text', y = .5, x = .25, label = '50%', size = 10) + 
  annotate('text', y = .95, x = .25, label = '100%', size = 10) + 
  annotate('text', y = .6, x = 1, size = 10,
           label = bquote(.(round(part.SUeff$part$fract[7,'Adj.R.square'], digits = 2) * 100) ~ "%"))

p.SUrange = ggplot(data = NULL) +
  geom_bar(aes(y = 1, x = 1), stat = 'identity', color = 'black', fill = 'white') +
  geom_bar(aes(y = part.SUrange$part$fract[3,'Adj.R.squared'], x = 1), stat = 'identity', color = 'black', fill = 'black') +
  coord_flip(ylim = c(0,1), xlim = c(.2,1.5), clip = 'off') +
  theme_classic() +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(), 
        panel.border = element_blank()) + 
  annotate('line', y = c(0,1), x = c(.45,.45)) +  
  annotate('line', y = c(0,0), x = c(.45,.4)) +
  annotate('line', y = c(.25,.25), x = c(.45,.4)) +  
  annotate('line', y = c(.5,.5), x = c(.45,.4)) +
  annotate('line', y = c(.75,.75), x = c(.45,.4)) +
  annotate('line', y = c(1,1), x = c(.45,.4)) +  
  annotate('text', y = 0.05, x = .25, label = '0%', size = 10) + 
  annotate('text', y = .5, x = .25, label = '50%', size = 10) + 
  annotate('text', y = .95, x = .25, label = '100%', size = 10) +
  annotate('text', y = .51, x = 1, size = 10,
           label = bquote(.(round(part.SUrange$part$fract[3,'Adj.R.squared'], digits = 2) * 100) ~ "%"))

ggsave(p.resp, filename = 'Fig3-var-resp.tiff', height = 6, width = 12.1, units = 'cm')
ggsave(p.SUrange, filename = 'Fig3-var-SUrange.tiff', height = 6, width = 12.1, units = 'cm')
ggsave(p.SUeff, filename = 'Fig3-var-SUeff.tiff', height = 6, width = 12.1, units = 'cm')

#### >> 2.3 Plot of microbial facets effects ####
df.plot$estimate[is.na(df.plot$estimate)] = 0
df.plot$exp = df.plot$exp %>% 
  factor(.,levels = c('FG.evenness','Cata',
                      'shannon.diversity.fungi', 'shannon.diversity.bacteria','bf',
                      'active.biomass','total.biomass'))

df.plot$resp = df.plot$resp %>% 
  factor(., levels = c('SU.eff', 'SU.range', 'microbial.respiration'))

df.plot$sig = NA
df.plot$sig[df.plot$pv<0.1& df.plot$pv<= 0.05] = '.'
df.plot$sig[df.plot$pv<0.05& df.plot$pv<= 0.01] = '*'
df.plot$sig[df.plot$pv<0.01& df.plot$pv<= 0.001] = '**'
df.plot$sig[df.plot$pv<0.001] = '***'

# plot
p.est = 
  ggplot(data = df.plot, 
         na.rm = TRUE, 
         aes(x = resp, y = exp , 
             color = estimate)) + 
  
  geom_point(data = df.plot, 
             aes(x = resp, y = exp,
                 color = estimate), 
             size = 20, pch = 15) + 
  
  annotate(geom = 'text', x = df.plot$resp, y = df.plot$exp, label = df.plot$sig, size = 10) + 
  
  scale_color_gradient2(low = 'Darkred', high = 'Darkblue', limit = c(-1, 1), guide = F) +
  theme(
    panel.grid = element_blank(),
    plot.background = element_rect(fill = 'white'),
    panel.background = element_rect(fill = 'white', colour = 'gray'),
    legend.position = 'right',
    plot.margin = margin(t = .8, r = .2, b = .6, l = 2, "cm")) + 
  theme_void() +
  coord_cartesian(ylim = c(1,8), xlim = c(1,3), clip = 'off')

ggsave(p.est, filename = 'Results/vf3/Figure3-drivers-1-3.tiff', height = 14, width = 15, units = 'cm')


#### Figure 3 ####
# The final figure combining all the plot have been made manually using Inkscape