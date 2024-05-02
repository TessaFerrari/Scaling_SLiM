# LOAD LIBRARIES AND READ IN DATA ----

# Load libraries
library(tidyverse)
library(scales)
library(ggpubr)
library(extrafont)
library(ggh4x)
#font_import(path="C:/Users/tferrari/AppData/Local/Microsoft/Windows/Fonts")
loadfonts(device = "win", quiet = TRUE)

# Set working directory
setwd("C:/Users/tferrari/Desktop/SlimBenchmark/fly/data")

# Read in data
mutStatData <- read_delim("mutStats_table.txt",delim="\t")


# CONVERT TICK OF ORIGIN INTO MUTATION AGE ----

# Subtract mutation tick from end-of-simulation tick
mutStatData$TotMutAge <- (mutStatData$CurrentTick - mutStatData$Tot.MeanTO)
mutStatData$DelMutAge <- (mutStatData$CurrentTick - mutStatData$Del.MeanTO)
mutStatData$NeuMutAge <- (mutStatData$CurrentTick - mutStatData$Neu.MeanTO)


# COMBINE EVERY 10 100kb GENOMES INTO 1 1Mb GENOME ----

mutStatData2 <- mutStatData 
mutStatData2$SimNum <- case_when(mutStatData2$GenomeSize==1e5 ~ floor((mutStatData2$SimNum-1)/10)+1,
                             mutStatData2$GenomeSize!=1e5 ~ mutStatData2$SimNum)
mutStatData2 <- mutStatData2 %>%
  group_by(BurnInType,ScalingFactor,GenomeSize,DomCoefficent,BurnNum,SimNum) %>%
  summarise(TotMeanS = mean(Tot.MeanS),
            DelMeanS = mean(Del.MeanS),
            TotMeanAge = mean(TotMutAge),
            DelMeanAge = mean(DelMutAge),
            NeuMeanAge = mean(NeuMutAge))


# SUMMARIZE REPLICATES WITH MEAN AND SD ----

mutStatSummary <- mutStatData2 %>%
  group_by(BurnInType,ScalingFactor,GenomeSize,DomCoefficent) %>%
  summarise(TotS.mean = mean(TotMeanS), TotS.sd = sd(TotMeanS), 
            DelS.mean = mean(DelMeanS), DelS.sd = sd(DelMeanS), 
            TotAge.mean = mean(TotMeanAge), TotAge.sd = sd(TotMeanAge),
            DelAge.mean = mean(DelMeanAge), DelAge.sd = sd(DelMeanAge),
            NeuAge.mean = mean(NeuMeanAge), NeuAge.sd = sd(NeuMeanAge))


# CALCULATE EXPECTED SELECTION COEFFICIENT ----

# Selection coefficient of mutations being added to the simulation (BEFORE mutations are filtered by selection)
mutStatSummary$ExpectedS <- (-0.000133*mutStatSummary$ScalingFactor)


# REFORMAT DATA ----

# Turn stats to numeric
mutStatSummary[,5:15] <- lapply(mutStatSummary[,5:15], as.numeric)

# Rename parameters and turn into factors
mutStatSummary$ScalingFactor <- as.factor(as.numeric(mutStatSummary$ScalingFactor))
mutStatSummary$BurnInType <- case_when(mutStatSummary$BurnInType=="5" ~ "5N",
                                   mutStatSummary$BurnInType=="10" ~ "10N",
                                   mutStatSummary$BurnInType=="20" ~ "20N",
                                   mutStatSummary$BurnInType=="Coal" ~ "Coal",
                                   mutStatSummary$BurnInType=="Recap" ~ "Recap")
mutStatSummary$BurnInType <- factor(mutStatSummary$BurnInType, levels=c('5N','10N','20N','Coal','Recap'))
mutStatSummary$GenomeSize <- case_when(mutStatSummary$GenomeSize==1e+05 ~ "100kb x 10",
                                   mutStatSummary$GenomeSize==1e+06 ~ "1Mb",
                                   mutStatSummary$GenomeSize==1e+07 ~ "10Mb")
mutStatSummary$GenomeSize <- factor(mutStatSummary$GenomeSize, levels=c('100kb x 10','1Mb','10Mb'))
mutStatSummary$DomCoefficent <- case_when(mutStatSummary$DomCoefficent==0.0 ~ "Recessive",
                                      mutStatSummary$DomCoefficent==0.5 ~ "Additive")
mutStatSummary$DomCoefficent <- factor(mutStatSummary$DomCoefficent, levels=c('Recessive','Additive'))


# GGPLOT THEME FORMATTING FOR GRAPHS ----

# Reference colors
ponyoPalette= c("#4D413F","#5A7080","#288B9A","#E75B64","#DE7862","#D8AF37","#E8C4A2","#F8E7D3")

# Custom theme function
theme_ponyo <- function(){ 
  graphPalette= c("100"="#288B9A","500"="#E75B64","1000"="#D8AF37")
  font = "Myriad Pro"   #assign font family up front
  dodge = position_dodge(width = 0.9)
  list(scale_fill_manual(name = "Scaling Factor", values = graphPalette),
       #scale_y_continuous(expand = expansion(mult = c(0, .05)),labels = function(x) format(x, scientific = TRUE)),
       geom_bar(stat = "identity", width=0.9, size = 0.5, position = dodge, show.legend = TRUE, colour="black"),
       theme_bw() %+replace%    #replace elements we want to change
         theme(
           #grid elements
           strip.background = element_rect(fill = "#F8E7D3", colour = "black"),
           #text elements
           axis.text.x = element_text(hjust = 0.5, vjust = 0.5, size = 14, family=font), 
           axis.text.y = element_text(size = 14, family=font), 
           plot.title = element_text(size = 18, hjust = 0.5, margin=margin(0,0,15,0), family=font), 
           axis.title = element_text(size = 16, family=font),
           strip.text = element_text(size = 14, family=font, margin=margin(5,5,5,5))
         ))
}
options(scipen=999) # Remove scientific notation


# SELECTION COEFFICIENT BAR PLOTS ----

# Total selection coefficient bar plot (including neutral)
TotS_afr <- ggplot(data = mutStatSummary, aes(x=BurnInType, y=TotS.mean, fill=ScalingFactor))  + 
  theme_ponyo() +
  geom_errorbar(aes(ymin=TotS.mean-TotS.sd, ymax=TotS.mean+TotS.sd), width=.2, position=position_dodge(.9)) +
  ggtitle("Mean Selection Coefficient of Fly African Population") +
  labs(x="Burn-In Method",y="Mean Selection Coefficient") + 
  facet_grid2(DomCoefficent ~ GenomeSize)+
  scale_y_reverse(limits = c(0,-0.1),expand = expansion(mult = c(0, .05)))

# Deleterious mutation selection coefficient bar plot
DelS_afr <- ggplot(data = mutStatSummary, aes(x=BurnInType, y=DelS.mean, fill=ScalingFactor))  + 
  theme_ponyo() +
  geom_errorbar(aes(ymin=DelS.mean-DelS.sd, ymax=DelS.mean+DelS.sd), width=.2, position=position_dodge(.9)) +
  geom_point(aes(x=BurnInType, y=ExpectedS), size = 2, position=position_dodge(.9), show.legend = FALSE) +
  ggtitle("Deleterious Mean Selection Coefficient of Fly African Population") +
  labs(x="Burn-In Method",y="Mean Selection Coefficient") + 
  facet_grid2(DomCoefficent ~ GenomeSize)+
  scale_y_reverse(limits = c(0,-0.135),expand = expansion(mult = c(0, .05)))


# MUTATION AGE BAR PLOTS ----

# Total mutation age bar plot (including neutral)
TotAge_afr <- ggplot(data = mutStatSummary, aes(x=BurnInType, y=TotAge.mean, fill=ScalingFactor))  + 
  theme_ponyo() +
  geom_errorbar(aes(ymin=TotAge.mean-TotAge.sd, ymax=TotAge.mean+TotAge.sd), width=.2, position=position_dodge(.9)) +
  ggtitle("Mean Mutation Age of Fly African Population") +
  labs(x="Burn-In Method",y="Mutation Age (ticks)") + 
  facet_grid2(DomCoefficent ~ GenomeSize, scales = "free_y", independent = "y")+
  scale_y_continuous(limits = c(0,500),expand = expansion(mult = c(0, .05)))

# Deleterious mutation age bar plot
DelAge_afr <- ggplot(data = mutStatSummary, aes(x=BurnInType, y=DelAge.mean, fill=ScalingFactor))  + 
  theme_ponyo() +
  geom_errorbar(aes(ymin=DelAge.mean-DelAge.sd, ymax=DelAge.mean+DelAge.sd), width=.2, position=position_dodge(.9)) +
  ggtitle("Mean Deleterious Mutation Age of Fly African Population") +
  labs(x="Burn-In Method",y="Mutation Age (ticks)") + 
  facet_grid2(DomCoefficent ~ GenomeSize)+
  scale_y_continuous(limits = c(0,275),expand = expansion(mult = c(0, .05)))

# Neutral mutation age bar plot (traditional burn ins)
NeuAge_burn_afr <- ggplot(data = mutStatSummary[(mutStatSummary$BurnInType!="Coal"&mutStatSummary$BurnInType!="Recap"),], aes(x=BurnInType, y=NeuAge.mean, fill=ScalingFactor))  + 
  theme_ponyo() +
  geom_errorbar(aes(ymin=NeuAge.mean-NeuAge.sd, ymax=NeuAge.mean+NeuAge.sd), width=.2, position=position_dodge(.9)) +
  ggtitle("Mean Neutral Mutation Age of Fly African Population") +
  labs(x="Burn-In Method",y="Mutation Age (ticks)") + 
  facet_grid2(DomCoefficent ~ GenomeSize)+
  scale_y_continuous(limits = c(0,530),expand = expansion(mult = c(0, .05)))

# Neutral mutation age bar plot (coalescence checking burn ins)
NeuAge_coal_afr <- ggplot(data = mutStatSummary[mutStatSummary$BurnInType=="Coal",], aes(x=BurnInType, y=NeuAge.mean, fill=ScalingFactor))  + 
  theme_ponyo() +
  geom_errorbar(aes(ymin=NeuAge.mean-NeuAge.sd, ymax=NeuAge.mean+NeuAge.sd), width=.2, position=position_dodge(.9)) +
  ggtitle("Mean Neutral Mutation Age of Fly African Population") +
  labs(x="Burn-In Method",y="Mutation Age (ticks)") + 
  facet_grid2(DomCoefficent ~ GenomeSize)+
  scale_y_continuous(limits = c(0,240000),expand = expansion(mult = c(0, .05)))

# Neutral mutation age bar plot (recapitated)
NeuAge_recap_afr <- ggplot(data = mutStatSummary[mutStatSummary$BurnInType=="Recap",], aes(x=BurnInType, y=NeuAge.mean, fill=ScalingFactor))  + 
  theme_ponyo() +
  geom_errorbar(aes(ymin=NeuAge.mean-NeuAge.sd, ymax=NeuAge.mean+NeuAge.sd), width=.2, position=position_dodge(.9)) +
  ggtitle("Mean Neutral Mutation Age of Fly African Population") +
  labs(x="Burn-In Method",y="Mutation Age (ticks)") + 
  facet_grid2(DomCoefficent ~ GenomeSize)+
  scale_y_continuous(limits = c(0,25000),expand = expansion(mult = c(0, .05)))


# SAVE PLOTS TO PDF ----

# Save selection coefficient plots
ggsave(TotS_afr, filename = "C:/Users/tferrari/Desktop/SlimBenchmark/figures/mutStatCounts/fly_TotS_graph.pdf", 
       device = cairo_pdf, width = 11, height = 8.5, units = "in")
ggsave(DelS_afr, filename = "C:/Users/tferrari/Desktop/SlimBenchmark/figures/mutStatCounts/fly_DelS_graph.pdf", 
       device = cairo_pdf, width = 11, height = 8.5, units = "in")

# Save mutation age plots
ggsave(TotAge_afr, filename = "C:/Users/tferrari/Desktop/SlimBenchmark/figures/mutStatCounts/fly_TotA_graph.pdf", 
       device = cairo_pdf, width = 11, height = 8.5, units = "in")
ggsave(TotDel_afr, filename = "C:/Users/tferrari/Desktop/SlimBenchmark/figures/mutStatCounts/fly_DelA_graph.pdf", 
       device = cairo_pdf, width = 11, height = 8.5, units = "in")
ggsave(TotNeu_burn_afr, filename = "C:/Users/tferrari/Desktop/SlimBenchmark/figures/mutStatCounts/fly_NeuBurnA_graph.pdf", 
       device = cairo_pdf, width = 11, height = 8.5, units = "in")
ggsave(TotNeu_coal_afr, filename = "C:/Users/tferrari/Desktop/SlimBenchmark/figures/mutStatCounts/fly_NeuCoalA_graph.pdf", 
       device = cairo_pdf, width = 11, height = 8.5, units = "in")
ggsave(TotNeu_recap_afr, filename = "C:/Users/tferrari/Desktop/SlimBenchmark/figures/mutStatCounts/fly_NeuRecapA_graph.pdf", 
       device = cairo_pdf, width = 11, height = 8.5, units = "in")
