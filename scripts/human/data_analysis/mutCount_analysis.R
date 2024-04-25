# Load libraries
library(tidyverse)
library(scales)
library(ggpubr)
library(ghibli)
library(extrafont)
library(ggh4x)
#font_import(path="C:/Users/tferrari/AppData/Local/Microsoft/Windows/Fonts")
loadfonts(device = "win", quiet = TRUE)

# Set working dir
setwd("C:/Users/tferrari/Desktop/SlimBenchmark/human/data")

# Read in data
mutData <- read_delim("mutCounts_table.txt",delim="\t")

# Combine every 10 100kb genomes into 1 1Mb genome
mutData2 <- mutData 
mutData2$SimNum <- case_when(mutData2$GenomeSize==1e5 ~ floor((mutData2$SimNum-1)/10)+1,
                             mutData2$GenomeSize!=1e5 ~ mutData2$SimNum)
mutData2 <- mutData2 %>%
  group_by(BurnInType,ScalingFactor,GenomeSize,DomCoefficent,BurnNum,SimNum) %>%
  summarise(ANeu = sum(ANeu), ENeu = sum(ENeu), EANeu = sum(EANeu), 
            ADel = sum(ADel), EDel = sum(EDel), EADel = sum(EADel))

# Add total mutation counts
mutData2$ATot <- mutData2$ADel + mutData2$ANeu
mutData2$ETot <- mutData2$EDel + mutData2$ENeu
mutData2$EATot <- mutData2$EADel + mutData2$EANeu

# Summarize replicates with mean and SD
mutSummary <- mutData2 %>%
  group_by(BurnInType,ScalingFactor,GenomeSize,DomCoefficent) %>%
  summarise(A.mean.neu = mean(ANeu), E.mean.neu = mean(ENeu), EA.mean.neu = mean(EANeu), 
            A.mean.del = mean(ADel), E.mean.del = mean(EDel), EA.mean.del = mean(EADel),
            A.mean.tot = mean(ATot), E.mean.tot = mean(ETot), EA.mean.tot = mean(EATot),
            A.sd.neu = sd(ANeu), E.sd.neu = sd(ENeu), EA.sd.neu = sd(EANeu), 
            A.sd.del = sd(ADel), E.sd.del = sd(EDel), EA.sd.del = sd(EADel),
            A.sd.tot = sd(ATot), E.sd.tot = sd(ETot), EA.sd.tot = sd(EATot))

# Reformat data
mutSummary[,5:22] <- lapply(mutSummary[,5:22], as.numeric)

mutSummary$ScalingFactor <- as.factor(as.numeric(mutSummary$ScalingFactor))

mutSummary$BurnInType <- case_when(mutSummary$BurnInType=="5" ~ "5N",
                                   mutSummary$BurnInType=="10" ~ "10N",
                                   mutSummary$BurnInType=="20" ~ "20N",
                                   mutSummary$BurnInType=="Coal" ~ "Coal",
                                   mutSummary$BurnInType=="Recap" ~ "Recap")
mutSummary$BurnInType <- factor(mutSummary$BurnInType, levels=c('5N','10N','20N','Coal','Recap'))

mutSummary$GenomeSize <- case_when(mutSummary$GenomeSize==1e+05 ~ "100kb x 10",
                                   mutSummary$GenomeSize==1e+06 ~ "1Mb",
                                   mutSummary$GenomeSize==1e+07 ~ "10Mb")
mutSummary$GenomeSize <- factor(mutSummary$GenomeSize, levels=c('100kb x 10','1Mb','10Mb'))

mutSummary$DomCoefficent <- case_when(mutSummary$DomCoefficent==0.0 ~ "Recessive",
                                      mutSummary$DomCoefficent==0.5 ~ "Additive")
mutSummary$DomCoefficent <- factor(mutSummary$DomCoefficent, levels=c('Recessive','Additive'))

# Graph formatting
ponyoPalette= c("#4D413F","#5A7080","#288B9A","#E75B64","#DE7862","#D8AF37","#E8C4A2","#F8E7D3")
theme_ponyo <- function(){ 
  graphPalette= c("1"="#288B9A","5"="#E75B64","10"="#D8AF37")
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

# African Population Neutral Mutation Counts Bar plot
neu_afr <- ggplot(data = mutSummary, aes(x=BurnInType, y=A.mean.neu, fill=ScalingFactor))  + 
  theme_ponyo() +
  geom_errorbar(aes(ymin=A.mean.neu-A.sd.neu, ymax=A.mean.neu+A.sd.neu), width=.2, position=position_dodge(.9)) +
  ggtitle("Neutral Mutation Counts of Human African Population") +
  labs(x="Burn-In Method",y="Neutral Mutation Counts") + 
  facet_grid2(DomCoefficent ~ GenomeSize, scales = "free_y", independent = "y")+
  scale_y_continuous(expand = expansion(mult = c(0, .05)))

# European Population Neutral Mutation Counts Bar Plot
neu_eur <- ggplot(data = mutSummary, aes(x=BurnInType, y=E.mean.neu, fill=ScalingFactor))  + 
  theme_ponyo() +
  geom_errorbar(aes(ymin=E.mean.neu-E.sd.neu, ymax=E.mean.neu+E.sd.neu), width=.2, position=position_dodge(.9)) +
  ggtitle("Neutral Mutation Counts of Human European Population") +
  labs(x="Burn-In Method",y="Neutral Mutation Counts") + 
  facet_grid2(DomCoefficent ~ GenomeSize, scales = "free_y", independent = "y")+
  scale_y_continuous(expand = expansion(mult = c(0, .05)))

# East Asian Population Neutral Mutation Counts Bar Plot
neu_easi <- ggplot(data = mutSummary, aes(x=BurnInType, y=EA.mean.neu, fill=ScalingFactor))  + 
  theme_ponyo() +
  geom_errorbar(aes(ymin=EA.mean.neu-EA.sd.neu, ymax=EA.mean.neu+EA.sd.neu), width=.2, position=position_dodge(.9)) +
  ggtitle("Neutral Mutation Counts of Human East Asian Population") +
  labs(x="Burn-In Method",y="Neutral Mutation Counts") + 
  facet_grid2(DomCoefficent ~ GenomeSize, scales = "free_y", independent = "y")+
  scale_y_continuous(expand = expansion(mult = c(0, .05)))

# African Population Deleterious Mutation Counts Bar plot
del_afr <- ggplot(data = mutSummary, aes(x=BurnInType, y=A.mean.del, fill=ScalingFactor))  + 
  theme_ponyo() +
  geom_errorbar(aes(ymin=A.mean.del-A.sd.del, ymax=A.mean.del+A.sd.del), width=.2, position=position_dodge(.9)) +
  ggtitle("Deleterious Mutation Counts of Human African Population") +
  labs(x="Burn-In Method",y="Deleterious Mutation Counts") + 
  facet_grid2(DomCoefficent ~ GenomeSize, scales = "free_y", independent = "y")+
  scale_y_continuous(expand = expansion(mult = c(0, .05)))

# European Population Deleterious Mutation Counts Bar plot
del_eur <- ggplot(data = mutSummary, aes(x=BurnInType, y=E.mean.del, fill=ScalingFactor))  + 
  theme_ponyo() +
  geom_errorbar(aes(ymin=E.mean.del-E.sd.del, ymax=E.mean.del+E.sd.del), width=.2, position=position_dodge(.9)) +
  ggtitle("Deleterious Mutation Counts of Human European Population") +
  labs(x="Burn-In Method",y="Deleterious Mutation Counts") + 
  facet_grid2(DomCoefficent ~ GenomeSize, scales = "free_y", independent = "y")+
  scale_y_continuous(expand = expansion(mult = c(0, .05)))

# East Asian Population Deleterious Mutation Counts Bar plot
del_easi <- ggplot(data = mutSummary, aes(x=BurnInType, y=EA.mean.del, fill=ScalingFactor))  + 
  theme_ponyo() +
  geom_errorbar(aes(ymin=EA.mean.del-EA.sd.del, ymax=EA.mean.del+EA.sd.del), width=.2, position=position_dodge(.9)) +
  ggtitle("Deleterious Mutation Counts of Human East Asian Population") +
  labs(x="Burn-In Method",y="Deleterious Mutation Counts") + 
  facet_grid2(DomCoefficent ~ GenomeSize, scales = "free_y", independent = "y")+
  scale_y_continuous(expand = expansion(mult = c(0, .05)))

# African Population Total Mutation Counts Bar plot
tot_afr <- ggplot(data = mutSummary, aes(x=BurnInType, y=A.mean.tot, fill=ScalingFactor))  + 
  theme_ponyo() +
  geom_errorbar(aes(ymin=A.mean.tot-A.sd.tot, ymax=A.mean.tot+A.sd.tot), width=.2, position=position_dodge(.9)) +
  ggtitle("Total Mutation Counts of Human African Population") +
  labs(x="Burn-In Method",y="Total Mutation Counts") + 
  facet_grid2(DomCoefficent ~ GenomeSize, scales = "free_y", independent = "y")+
  scale_y_continuous(expand = expansion(mult = c(0, .05)))

# European Population Total Mutation Counts Bar Plot
tot_eur <- ggplot(data = mutSummary, aes(x=BurnInType, y=E.mean.tot, fill=ScalingFactor))  + 
  theme_ponyo() +
  geom_errorbar(aes(ymin=E.mean.tot-E.sd.tot, ymax=E.mean.tot+E.sd.tot), width=.2, position=position_dodge(.9)) +
  ggtitle("Total Mutation Counts of Human European Population") +
  labs(x="Burn-In Method",y="Total Mutation Counts") + 
  facet_grid2(DomCoefficent ~ GenomeSize, scales = "free_y", independent = "y")+
  scale_y_continuous(expand = expansion(mult = c(0, .05)))

# East Asian Population Total Mutation Counts Bar Plot
tot_easi <- ggplot(data = mutSummary, aes(x=BurnInType, y=EA.mean.tot, fill=ScalingFactor))  + 
  theme_ponyo() +
  geom_errorbar(aes(ymin=EA.mean.tot-EA.sd.tot, ymax=EA.mean.tot+EA.sd.tot), width=.2, position=position_dodge(.9)) +
  ggtitle("Total Mutation Counts of Human East Asian Population") +
  labs(x="Burn-In Method",y="Total Mutation Counts") + 
  facet_grid2(DomCoefficent ~ GenomeSize, scales = "free_y", independent = "y")+
  scale_y_continuous(expand = expansion(mult = c(0, .05)))


# Save plots to pdfs
ggsave(neu_afr, filename = "C:/Users/tferrari/Desktop/SlimBenchmark/figures/human_AfrNeu_graph.pdf", 
       device = cairo_pdf, width = 11, height = 8.5, units = "in")
ggsave(neu_eur, filename = "C:/Users/tferrari/Desktop/SlimBenchmark/figures/human_EurNeu_graph.pdf", 
       device = cairo_pdf, width = 11, height = 8.5, units = "in")
ggsave(neu_easi, filename = "C:/Users/tferrari/Desktop/SlimBenchmark/figures/human_EAsiNeu_graph.pdf", 
       device = cairo_pdf, width = 11, height = 8.5, units = "in")
ggsave(del_afr, filename = "C:/Users/tferrari/Desktop/SlimBenchmark/figures/human_AfrDel_graph.pdf", 
       device = cairo_pdf, width = 11, height = 8.5, units = "in")
ggsave(del_eur, filename = "C:/Users/tferrari/Desktop/SlimBenchmark/figures/human_EurDel_graph.pdf", 
       device = cairo_pdf, width = 11, height = 8.5, units = "in")
ggsave(del_easi, filename = "C:/Users/tferrari/Desktop/SlimBenchmark/figures/human_EAsiDel_graph.pdf", 
       device = cairo_pdf, width = 11, height = 8.5, units = "in")
ggsave(tot_afr, filename = "C:/Users/tferrari/Desktop/SlimBenchmark/figures/human_AfrTot_graph.pdf", 
       device = cairo_pdf, width = 11, height = 8.5, units = "in")
ggsave(tot_eur, filename = "C:/Users/tferrari/Desktop/SlimBenchmark/figures/human_EurTot_graph.pdf", 
       device = cairo_pdf, width = 11, height = 8.5, units = "in")
ggsave(tot_easi, filename = "C:/Users/tferrari/Desktop/SlimBenchmark/figures/human_EAsiTot_graph.pdf", 
       device = cairo_pdf, width = 11, height = 8.5, units = "in")
