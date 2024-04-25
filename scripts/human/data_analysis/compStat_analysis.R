# Load libraries ----
library(tidyverse)
library(scales)
library(ggpubr)
library(ghibli)
library(extrafont)
#font_import(path="C:/Users/tferrari/AppData/Local/Microsoft/Windows/Fonts")
loadfonts(device = "win", quiet = TRUE)

# Set working dir ----
setwd("C:/Users/tferrari/Desktop/SlimBenchmark/human/data")

# Read in data ----
burn_compData <- read_delim("burnin_compStat_table.tsv",delim="\t") 
main_compData <- read_delim("main_compStat_table.txt",delim="\t")
msprime_compData <- read_delim("msprime_compStat_table.txt",delim="\t") 

# Combine every 10 100kb genomes into 1 1Mb genome ----
main_compData2 <- main_compData 
main_compData2$SimNum <- case_when(main_compData2$GenomeSize==1e5 ~ floor((main_compData2$SimNum-1)/10)+1,
                                   main_compData2$GenomeSize!=1e5 ~ main_compData2$SimNum)
main_compData2 <- main_compData2 %>%
  group_by(BurnInType,ScalingFactor,GenomeSize,DomCoefficent,BurnNum,SimNum) %>%
  summarise(MemMB=max(MemMB), TimeSecs=max(TimeSecs))

msprime_compData2 <- msprime_compData 
msprime_compData2$SimNum <- case_when(msprime_compData2$GenomeSize==1e5 ~ floor((msprime_compData2$SimNum-1)/10)+1,
                                      msprime_compData2$GenomeSize!=1e5 ~ msprime_compData2$SimNum)
msprime_compData2 <- msprime_compData2 %>%
  group_by(BurnInType,ScalingFactor,GenomeSize,DomCoefficent,BurnNum,SimNum) %>%
  summarise(TimeSecs=max(TimeSecs))


# Combine data (come back to this later, need to figure out how to incorporate ----
#               msprime runtime into the total simulation runtimes)
tot_compData <- merge(burn_compData, main_compData2, all=TRUE,
                      by.x=c("BurnInType","ScalingFactor","GenomeSize","RepNum"),
                      by.y=c("BurnInType","ScalingFactor","GenomeSize","BurnNum"),
                      suffixes = c("_burn","_main"))
tot_compData <- merge(tot_compData, msprime_compData2, all=TRUE,
                      by.x=c("BurnInType","ScalingFactor","GenomeSize","DomCoefficent","RepNum","SimNum"),
                      by.y=c("BurnInType","ScalingFactor","GenomeSize","DomCoefficent","BurnNum","SimNum"),
                      suffixes = c("","_msprime"))
tot_compData <- tot_compData %>% rename(TimeSecs_msprime=TimeSecs)
tot_compData <- replace(tot_compData, is.na(tot_compData), 0)
tot_compData$TimeSecs_tot <- tot_compData$TimeSecs_burn+tot_compData$TimeSecs_main+tot_compData$TimeSecs_msprime
tot_compData <- tot_compData[,c(1:6,8,10:12)]

# Create tree generation data ----
tree_compData <- tot_compData %>%
  group_by(BurnInType,ScalingFactor,GenomeSize,RepNum) %>%
  summarise(TimeSecs_tree = mean(TimeSecs_burn))
tree_compData <- tree_compData[tree_compData$BurnInType!="Coal" & tree_compData$BurnInType!="Recap",]
temp <- tot_compData[tot_compData$BurnInType=="Coal" | tot_compData$BurnInType=="Recap",c(1:7,9)]
temp$TimeSecs_tree <- temp$TimeSecs_burn + temp$TimeSecs_msprime
tree_compData <- merge(tree_compData, temp, all =TRUE,
                      by.x=c("BurnInType","ScalingFactor","GenomeSize","RepNum","TimeSecs_tree"),
                      by.y=c("BurnInType","ScalingFactor","GenomeSize","RepNum","TimeSecs_tree"),
                      suffixes = c("","_msprime"))[,1:7]


# Summarize replicates with mean and SD ----
burn_compSummary <- burn_compData %>%
  group_by(BurnInType,ScalingFactor,GenomeSize) %>%
  summarise(meanRT = mean(TimeSecs), sdRT = sd(TimeSecs),
            meanMU = mean(MemMB), sdMU = sd(MemMB))
main_compSummary <- main_compData2 %>%
  group_by(BurnInType,ScalingFactor,GenomeSize,DomCoefficent) %>%
  summarise(meanRT = mean(TimeSecs), sdRT = sd(TimeSecs),
            meanMU = mean(MemMB), sdMU = sd(MemMB))
tot_compSummary <- tot_compData %>%
  group_by(BurnInType,ScalingFactor,GenomeSize,DomCoefficent) %>%
  summarise(meanRT = mean(TimeSecs_tot), sdRT = sd(TimeSecs_tot))
tree_compSummary <- tree_compData %>%
  group_by(BurnInType,ScalingFactor,GenomeSize) %>%
  summarise(meanRT = mean(TimeSecs_tree), sdRT = sd(TimeSecs_tree))

# Reformat burn in data ----
burn_compSummary[,4:7] <- lapply(burn_compSummary[,4:7], as.numeric)
burn_compSummary$ScalingFactor <- as.factor(as.numeric(burn_compSummary$ScalingFactor))
burn_compSummary$BurnInType <- case_when(burn_compSummary$BurnInType=="5" ~ "5N",
                                   burn_compSummary$BurnInType=="10" ~ "10N",
                                   burn_compSummary$BurnInType=="20" ~ "20N",
                                   burn_compSummary$BurnInType=="Coal" ~ "Coal",
                                   burn_compSummary$BurnInType=="Recap" ~ "Recap")
burn_compSummary$BurnInType <- factor(burn_compSummary$BurnInType, levels=c('5N','10N','20N','Coal','Recap'))
burn_compSummary$GenomeSize <- case_when(burn_compSummary$GenomeSize==1e+05 ~ "100kb x 10",
                                   burn_compSummary$GenomeSize==1e+06 ~ "1Mb",
                                   burn_compSummary$GenomeSize==1e+07 ~ "10Mb")
burn_compSummary$GenomeSize <- factor(burn_compSummary$GenomeSize, levels=c('100kb x 10','1Mb','10Mb'))

# Reformat main data ----
main_compSummary[,5:8] <- lapply(main_compSummary[,5:8], as.numeric)
main_compSummary$ScalingFactor <- as.factor(as.numeric(main_compSummary$ScalingFactor))
main_compSummary$BurnInType <- case_when(main_compSummary$BurnInType=="5" ~ "5N",
                                         main_compSummary$BurnInType=="10" ~ "10N",
                                         main_compSummary$BurnInType=="20" ~ "20N",
                                         main_compSummary$BurnInType=="Coal" ~ "Coal",
                                         main_compSummary$BurnInType=="Recap" ~ "Recap")
main_compSummary$BurnInType <- factor(main_compSummary$BurnInType, levels=c('5N','10N','20N','Coal','Recap'))
main_compSummary$GenomeSize <- case_when(main_compSummary$GenomeSize==1e+05 ~ "100kb x 10",
                                         main_compSummary$GenomeSize==1e+06 ~ "1Mb",
                                         main_compSummary$GenomeSize==1e+07 ~ "10Mb")
main_compSummary$GenomeSize <- factor(main_compSummary$GenomeSize, levels=c('100kb x 10','1Mb','10Mb'))
main_compSummary$DomCoefficent <- case_when(main_compSummary$DomCoefficent==0.0 ~ "Recessive",
                                      main_compSummary$DomCoefficent==0.5 ~ "Additive")
main_compSummary$DomCoefficent <- factor(main_compSummary$DomCoefficent, levels=c('Recessive','Additive'))

# Reformat total data ----
tot_compSummary[,5:6] <- lapply(tot_compSummary[,5:6], as.numeric)
tot_compSummary$ScalingFactor <- as.factor(as.numeric(tot_compSummary$ScalingFactor))
tot_compSummary$BurnInType <- case_when(tot_compSummary$BurnInType=="5" ~ "5N",
                                         tot_compSummary$BurnInType=="10" ~ "10N",
                                         tot_compSummary$BurnInType=="20" ~ "20N",
                                         tot_compSummary$BurnInType=="Coal" ~ "Coal",
                                         tot_compSummary$BurnInType=="Recap" ~ "Recap")
tot_compSummary$BurnInType <- factor(tot_compSummary$BurnInType, levels=c('5N','10N','20N','Coal','Recap'))
tot_compSummary$GenomeSize <- case_when(tot_compSummary$GenomeSize==1e+05 ~ "100kb x 10",
                                         tot_compSummary$GenomeSize==1e+06 ~ "1Mb",
                                         tot_compSummary$GenomeSize==1e+07 ~ "10Mb")
tot_compSummary$GenomeSize <- factor(tot_compSummary$GenomeSize, levels=c('100kb x 10','1Mb','10Mb'))
tot_compSummary$DomCoefficent <- case_when(tot_compSummary$DomCoefficent==0.0 ~ "Recessive",
                                            tot_compSummary$DomCoefficent==0.5 ~ "Additive")
tot_compSummary$DomCoefficent <- factor(tot_compSummary$DomCoefficent, levels=c('Recessive','Additive'))

# Reformat tree data ----
tree_compSummary[,4:5] <- lapply(tree_compSummary[,4:5], as.numeric)
tree_compSummary$ScalingFactor <- as.factor(as.numeric(tree_compSummary$ScalingFactor))
tree_compSummary$BurnInType <- case_when(tree_compSummary$BurnInType=="5" ~ "5N",
                                        tree_compSummary$BurnInType=="10" ~ "10N",
                                        tree_compSummary$BurnInType=="20" ~ "20N",
                                        tree_compSummary$BurnInType=="Coal" ~ "Coal",
                                        tree_compSummary$BurnInType=="Recap" ~ "Recap")
tree_compSummary$BurnInType <- factor(tree_compSummary$BurnInType, levels=c('5N','10N','20N','Coal','Recap'))
tree_compSummary$GenomeSize <- case_when(tree_compSummary$GenomeSize==1e+05 ~ "100kb x 10",
                                        tree_compSummary$GenomeSize==1e+06 ~ "1Mb",
                                        tree_compSummary$GenomeSize==1e+07 ~ "10Mb")
tree_compSummary$GenomeSize <- factor(tree_compSummary$GenomeSize, levels=c('100kb x 10','1Mb','10Mb'))

# Graph formatting ----
ponyoPalette= c("#4D413F","#5A7080","#288B9A","#E75B64","#DE7862","#D8AF37","#E8C4A2","#F8E7D3")
theme_ponyo <- function(){ 
  graphPalette= c('1'="#288B9A",'5'="#E75B64",'10'="#D8AF37")
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

# Runtime bar plots ----
# Runtime Bar plot for burn in
#burn_runtime <- ggplot(data = burn_compSummary, aes(x=BurnInType, y=meanRT, fill=ScalingFactor))  + 
#  theme_ponyo() +
#  geom_errorbar(aes(ymin=meanRT-sdRT, ymax=meanRT+sdRT), width=.2, position=position_dodge(.9)) +  
#  ggtitle("Runtime of Human Burn-In Simulation") +
#  labs(x="Coalescence Method",y="Runtime (minutes)") + 
#  facet_wrap(. ~ GenomeSize) +
#  scale_y_log10(expand = expansion(mult = c(0, .05)),labels = label_number(suffix = "", scale = 1/60))

# Runtime plot for main
main_runtime <- ggplot(data = main_compSummary, aes(x=BurnInType, y=meanRT, fill=ScalingFactor))  + 
  theme_ponyo() +
  geom_errorbar(aes(ymin=meanRT-sdRT, ymax=meanRT+sdRT), width=.2, position=position_dodge(.9)) +  
  ggtitle("Runtime of Human Main Simulations") +
  labs(x="Coalescence Method",y="Runtime (minutes)") + 
  facet_grid(DomCoefficent ~ GenomeSize) +
  scale_y_log10(breaks=c(30,300,3000,30000),expand = expansion(mult = c(0, .05)),labels = label_number(suffix = "", scale = 1/60))

# Runtime plot for total runtime (burnin + main + msprime)
tot_runtime <- ggplot(data = tot_compSummary, aes(x=BurnInType, y=meanRT, fill=ScalingFactor))  + 
  theme_ponyo() +
  geom_errorbar(aes(ymin=meanRT-sdRT, ymax=meanRT+sdRT), width=.2, position=position_dodge(.9)) +  
  ggtitle("Total Runtime of Human Simulations") +
  labs(x="Coalescence Method",y="Runtime (minutes)") + 
  facet_grid(DomCoefficent ~ GenomeSize) +
  scale_y_log10(breaks=c(30,300,3000,30000),expand = expansion(mult = c(0, .05)),labels = label_number(suffix = "", scale = 1/60))

# Runtime plot for tree creation runtime (burnin + msprime)
tree_runtime <- ggplot(data = tree_compSummary, aes(x=BurnInType, y=meanRT, fill=ScalingFactor))  + 
  theme_ponyo() +
  geom_errorbar(aes(ymin=meanRT-sdRT, ymax=meanRT+sdRT), width=.2, position=position_dodge(.9)) +  
  ggtitle("Runtime of Human Simulation Coalescence") +
  labs(x="Coalescence Method",y="Runtime (minutes)") + 
  facet_wrap(. ~ GenomeSize) +
  scale_y_log10(breaks=c(0.3,3,30,300,3000),expand = expansion(mult = c(0.05, .05)),labels = label_number(suffix = "", scale = 1/60))

# Memory usage plots ----
# Memory Bar plot for burn in
burn_memory <- ggplot(data = burn_compSummary, aes(x=BurnInType, y=meanMU, fill=ScalingFactor))  + 
  theme_ponyo() +
  geom_errorbar(aes(ymin=meanMU-sdMU, ymax=meanMU+sdMU), width=.2, position=position_dodge(.9)) +  
  ggtitle("Memory Usage of Human Simulation Coalescence") +
  labs(x="Coalescence Method",y="Peak Memory Usage (MB)") + 
  facet_wrap(. ~ GenomeSize) +
  scale_y_log10(expand = expansion(mult = c(0, .05)))

# Memory plot for main
main_memory <- ggplot(data = main_compSummary, aes(x=BurnInType, y=meanMU, fill=ScalingFactor))  + 
  theme_ponyo() +
  geom_errorbar(aes(ymin=meanMU-sdMU, ymax=meanMU+sdMU), width=.2, position=position_dodge(.9)) +  
  ggtitle("Memory Usage of Human Main Simulations") +
  labs(x="Coalescence Method",y="Peak Memory Usage (MB)") + 
  facet_grid(DomCoefficent ~ GenomeSize) +
  scale_y_log10(expand = expansion(mult = c(0, .05)))



# Save plots to pdfs ----
#ggsave(burn_runtime, filename = "C:/Users/tferrari/Desktop/SlimBenchmark/figures/statsComp/human_burnRT_graph.pdf", 
#       device = cairo_pdf, width = 11, height = 8.5, units = "in")
ggsave(main_runtime, filename = "C:/Users/tferrari/Desktop/SlimBenchmark/figures/statsComp/human_mainRT_graph.pdf", 
       device = cairo_pdf, width = 11, height = 8.5, units = "in")
ggsave(tot_runtime, filename = "C:/Users/tferrari/Desktop/SlimBenchmark/figures/statsComp/human_totRT_graph.pdf", 
       device = cairo_pdf, width = 11, height = 8.5, units = "in")
ggsave(tree_runtime, filename = "C:/Users/tferrari/Desktop/SlimBenchmark/figures/statsComp/human_treeRT_graph.pdf", 
       device = cairo_pdf, width = 11, height = 8.5, units = "in")
ggsave(burn_memory, filename = "C:/Users/tferrari/Desktop/SlimBenchmark/figures/statsComp/human_burnMem_graph.pdf", 
       device = cairo_pdf, width = 11, height = 8.5, units = "in")
ggsave(main_memory, filename = "C:/Users/tferrari/Desktop/SlimBenchmark/figures/statsComp/human_mainMem_graph.pdf", 
       device = cairo_pdf, width = 11, height = 8.5, units = "in")
