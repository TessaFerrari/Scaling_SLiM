# Load libraries and read in data ----
library(tidyverse)
library(extrafont)
library(scales)
library(ggpubr)
library(ggh4x)
library(ggbreak) 
loadfonts(device = "win", quiet = TRUE)


setwd("C:/Users/tferrari/Desktop/SlimBenchmark/human/data")
sfsData <- as.data.frame(read_delim("subsamp50_sfs_table.txt",delim="\t"))

# Read in and reformat empirical SFS ----
# african data
multiIdx <- expand.grid(c("african"),c(5,10,20,"Coal","Recap"),c("Empirical"),c(1e5,1e6,1e7),c(0.0,0.5),c(1),c(1))
#YRIdata <- as.data.frame(t(read_delim("empiricalData/YRI_WholeGenomeCountsforSFS_Folded.txt",delim="\t")))[2,] # Whole genome data (totLength=2897310462)
YRIdata <- as.data.frame(t(read_delim("empiricalData/YRI_Exonic_FoldedN30.txt",delim="\t")))[2,] # Exonic data (totLength=26824075)
YRIdata <- rbind(YRIdata, YRIdata[rep(1, 29), ])
YRIdata <- cbind(multiIdx,YRIdata)


# european data
multiIdx <- expand.grid(c("european"),c(5,10,20,"Coal","Recap"),c("Empirical"),c(1e5,1e6,1e7),c(0.0,0.5),c(1),c(1))
#CEUdata <- as.data.frame(t(read_delim("empiricalData/CEU_WholeGenomeCountsforSFS_Folded.txt",delim="\t")))[2,] # Whole genome data (totLength=2897310462)
CEUdata <- as.data.frame(t(read_delim("empiricalData/CEU_Exonic_FoldedN30.txt",delim="\t")))[2,] # Exonic data (totLength=26824075)
CEUdata <- rbind(CEUdata, CEUdata[rep(1, 29), ])
CEUdata <- cbind(multiIdx,CEUdata)

# Combine YRI and CEU
empiricalData <- rbind(YRIdata,CEUdata)

# Reformat and add number of SNPs and callable sites
colnames(empiricalData) <- c("Population", "BurnInType", "ScalingFactor", "GenomeSize", "DomCoefficent", "BurnNum", "SimNum", as.character(0:30))
rownames(empiricalData) <-c(1:60)
empiricalData <- empiricalData %>% 
  rowwise() %>%
  mutate(TotNumMuts = sum(c_across(`1`:`30`)),
         TotSites = sum(c_across(`0`:`30`)),
         .after = 7) %>%
  select(-`0`)

# Sum every 10 100kb replicates to create 5 reps of 100kb x 10 ----
sfsData <- rbind.data.frame(sfsData[sfsData$GenomeSize==1e5,] %>% 
                              mutate(SimNum=floor((sfsData[sfsData$GenomeSize==1e5,]$SimNum-1)/10)+1) %>%
                              group_by(Population,BurnInType,ScalingFactor,GenomeSize,DomCoefficent,BurnNum,SimNum) %>%
                              summarise_at(as.character(c("TotNumMuts",1:50)), list(sum=sum)) %>%
                              rename_with(~ c("TotNumMuts",1:50), 8:58),
                            sfsData[sfsData$GenomeSize!=1e5,])

# Calculate population stats from simulation SFS ----

# Takes in freqs vector of segregating sites, returns heterozygosity (2pq) UNAJUSTED BY GENSIZE
calc_het <- function(freqs, nSamp) {
  pq = (1:nSamp/(2*nSamp))*(1-1:nSamp/(2*nSamp))
  exp_het <- sum(2*pq*freqs)
  return(exp_het)
}
# Takes in tot# of segregating sites, returns theta (K/a_n)) UNAJUSTED BY GENSIZE
calc_theta <-function(K,nSamp){
  a_n = sum(1/seq(1,(nSamp*2)-1,1))
  theta = K/a_n
  return(theta)
}
# Add popstats to dataframe
popsfsData <- sfsData %>% 
  rowwise() %>%
  mutate( ExpHet = calc_het(c_across(`1`:`50`),50)/GenomeSize, 
          Theta = calc_theta(TotNumMuts,50)/GenomeSize) %>%
  select(Population:SimNum, ExpHet, Theta)
# since 100kb are aggregated in sets of 10, divide by 10
popsfsData[popsfsData$GenomeSize==1e5,c("ExpHet","Theta")] = popsfsData[popsfsData$GenomeSize==1e5,c("ExpHet","Theta")]/10

# Calculate population stats from empirical SFS ----

# Add popstats to dataframe
popEmpiricalData <- empiricalData %>% 
  rowwise() %>%
  mutate( ExpHet = calc_het(c_across(`1`:`30`),30)/26824075, # divided by length of empirical data
          Theta = calc_theta(TotNumMuts,30)/26824075) %>%    # exonic=26824075, whole genome=2897310462
  select(Population:SimNum, ExpHet, Theta)

# add empirical data to pop stat data frame
popsfsData <- rbind(popsfsData,popEmpiricalData)

# Put Minor allele freq into bins (i.e. cols 1:5 sum to be MAF 1-5%) ----
sfsData$`0-0.05` <- rowSums(sfsData[,as.character(1:5)])
sfsData$`0.05-0.1` <- rowSums(sfsData[,as.character(6:10)])
sfsData$`0.1-0.15` <- rowSums(sfsData[,as.character(11:15)])
sfsData$`0.15-0.2` <- rowSums(sfsData[,as.character(16:20)])
sfsData$`0.2-0.25` <- rowSums(sfsData[,as.character(21:25)])
sfsData$`0.25-0.3` <- rowSums(sfsData[,as.character(26:30)])
sfsData$`0.3-0.35` <- rowSums(sfsData[,as.character(31:35)])
sfsData$`0.35-0.4` <- rowSums(sfsData[,as.character(36:40)])
sfsData$`0.4-0.45` <- rowSums(sfsData[,as.character(41:45)])
sfsData$`0.45-0.5` <- rowSums(sfsData[,as.character(46:50)])

# Bin by MAF for empirical data
empiricalData$`0-0.05` <- rowSums(empiricalData[,as.character(1:3)])
empiricalData$`0.05-0.1` <- rowSums(empiricalData[,as.character(4:6)])
empiricalData$`0.1-0.15` <- rowSums(empiricalData[,as.character(7:9)])
empiricalData$`0.15-0.2` <- rowSums(empiricalData[,as.character(10:12)])
empiricalData$`0.2-0.25` <- rowSums(empiricalData[,as.character(13:15)])
empiricalData$`0.25-0.3` <- rowSums(empiricalData[,as.character(16:18)])
empiricalData$`0.3-0.35` <- rowSums(empiricalData[,as.character(19:21)])
empiricalData$`0.35-0.4` <- rowSums(empiricalData[,as.character(22:24)])
empiricalData$`0.4-0.45` <- rowSums(empiricalData[,as.character(25:27)])
empiricalData$`0.45-0.5` <- rowSums(empiricalData[,as.character(28:30)])

# Remove small bins and combine simulated and empirical data---
sfsData <- sfsData[,c(1:8,59:68)]
empiricalData <- empiricalData[,c(1:8,40:49)]

sfsData$ScalingFactor <- as.character(sfsData$ScalingFactor)
sfsData <- rbind(sfsData,empiricalData)

# Make allele frequencies a proportion of total snp counts ----
sfsData[,9:18] <- sfsData[,9:18]/sfsData$TotNumMuts

# Add column for the sum of the sfs tail (cols 5 to 50)
#sfsData <- cbind.data.frame(sfsData[,1:12], sfsData[,13:58] %>% mutate(TailSum = rowSums(.)))

# Summarize replicates with mean and SD ----
sfsSummary <- sfsData %>%
  group_by(Population,BurnInType,ScalingFactor,GenomeSize,DomCoefficent) %>%
  summarise_at(as.character(c("0-0.05","0.05-0.1","0.1-0.15","0.15-0.2","0.2-0.25","0.25-0.3","0.3-0.35","0.35-0.4","0.4-0.45","0.45-0.5")), list(mean = mean, sd = sd)) %>% 
  replace(is.na(.), 0)

popsfsSummary <- popsfsData %>%
  group_by(Population,BurnInType,ScalingFactor,GenomeSize,DomCoefficent) %>%
  summarise_at(c("ExpHet","Theta"), list(mean = mean, sd = sd))  %>% 
  replace(is.na(.), 0)

# Reformat SFS Data ----
sfsSummary$ScalingFactor <- factor(sfsSummary$ScalingFactor, levels = c("Empirical","1","5","10"))
sfsSummary$Population <- case_when(sfsSummary$Population=="african" ~ "African",
                                   sfsSummary$Population=="european" ~ "European",
                                   sfsSummary$Population=="eastasian" ~ "East Asian")
sfsSummary$BurnInType <- case_when(sfsSummary$BurnInType=="5" ~ "5N",
                                   sfsSummary$BurnInType=="10" ~ "10N",
                                   sfsSummary$BurnInType=="20" ~ "20N",
                                   sfsSummary$BurnInType=="Coal" ~ "Coal",
                                   sfsSummary$BurnInType=="Recap" ~ "Recap")
sfsSummary$BurnInType <- factor(sfsSummary$BurnInType, levels=c('5N','10N','20N','Coal','Recap'))

sfsSummary$GenomeSize <- case_when(sfsSummary$GenomeSize==1e+05 ~ "100kb x 10",
                                   sfsSummary$GenomeSize==1e+06 ~ "1Mb",
                                   sfsSummary$GenomeSize==1e+07 ~ "10Mb")
sfsSummary$GenomeSize <- factor(sfsSummary$GenomeSize, levels=c('100kb x 10','1Mb','10Mb'))

sfsSummary$DomCoefficent <- case_when(sfsSummary$DomCoefficent==0.0 ~ "Recessive",
                                      sfsSummary$DomCoefficent==0.5 ~ "Additive")
sfsSummary$DomCoefficent <- factor(sfsSummary$DomCoefficent, levels=c('Recessive','Additive'))

# Reformat pop stat Data ----
popsfsSummary$ScalingFactor <- factor(popsfsSummary$ScalingFactor, levels = c("Empirical","1","5","10"))
popsfsSummary$Population <- case_when(popsfsSummary$Population=="african" ~ "African",
                                   popsfsSummary$Population=="european" ~ "European",
                                   popsfsSummary$Population=="eastasian" ~ "East Asian")
popsfsSummary$BurnInType <- case_when(popsfsSummary$BurnInType=="5" ~ "5N",
                                   popsfsSummary$BurnInType=="10" ~ "10N",
                                   popsfsSummary$BurnInType=="20" ~ "20N",
                                   popsfsSummary$BurnInType=="Coal" ~ "Coal",
                                   popsfsSummary$BurnInType=="Recap" ~ "Recap")
popsfsSummary$BurnInType <- factor(popsfsSummary$BurnInType, levels=c('5N','10N','20N','Coal','Recap'))

popsfsSummary$GenomeSize <- case_when(popsfsSummary$GenomeSize==1e+05 ~ "100kb x 10",
                                   popsfsSummary$GenomeSize==1e+06 ~ "1Mb",
                                   popsfsSummary$GenomeSize==1e+07 ~ "10Mb")
popsfsSummary$GenomeSize <- factor(popsfsSummary$GenomeSize, levels=c('100kb x 10','1Mb','10Mb'))

popsfsSummary$DomCoefficent <- case_when(popsfsSummary$DomCoefficent==0.0 ~ "Recessive",
                                      popsfsSummary$DomCoefficent==0.5 ~ "Additive")
popsfsSummary$DomCoefficent <- factor(popsfsSummary$DomCoefficent, levels=c('Recessive','Additive'))


# Graph formatting ----
ponyoPalette= c("#4D413F","#5A7080","#288B9A","#E75B64","#DE7862","#D8AF37","#E8C4A2","#F8E7D3")
theme_ponyo <- function(){ 
  graphPalette= c("Empirical"="#5A7080","1"="#288B9A","5"="#E75B64","10"="#D8AF37")
  font = "Myriad Pro"   #assign font family up front
  dodge = position_dodge(width = 0.9)
  list(scale_fill_manual(name = "Scaling Factor", values = graphPalette),
       #scale_y_continuous(expand = expansion(mult = c(0, .05)),labels = function(x) format(x, scientific = TRUE)),
       geom_bar(stat = "identity", width=0.9, size = 0.2, position = dodge, show.legend = TRUE, colour="black"),
       theme_bw() %+replace%    #replace elements we want to change
         theme(
           #grid elements
           strip.background = element_rect(fill = "#F8E7D3", colour = "black"),
           #text elements
           axis.text.x = element_text(hjust = 0.5, vjust = 0.5, size = 8, family=font), 
           axis.text.y = element_text(size = 12, family=font), 
           plot.title = element_text(size = 18, hjust = 0.5, margin=margin(0,0,15,0), family=font), 
           axis.title = element_text(size = 16, family=font),
           strip.text = element_text(size = 14, family=font, margin=margin(5,5,5,5))
         ))
}
options(scipen=10000) # remove scientific notation


# Pivot mean and sd ----
sfsLonger <- sfsSummary %>% 
  pivot_longer(cols = ends_with(c("_mean","_sd")), 
               names_to = c("AlleleFreq", ".value"), 
               names_sep="_" )%>% 
  filter(grepl( '-', AlleleFreq, fixed = TRUE))
sfsLonger$AlleleFreq <- factor(sfsLonger$AlleleFreq, levels = c("0-0.05","0.05-0.1","0.1-0.15","0.15-0.2","0.2-0.25","0.25-0.3","0.3-0.35","0.35-0.4","0.4-0.45","0.45-0.5"))

# Separate by population ----
sfs_afr <- sfsLonger[sfsLonger$Population=="African",]
sfs_eur <- sfsLonger[sfsLonger$Population=="European",]
sfs_easi <- sfsLonger[sfsLonger$Population=="East Asian",]

pop_afr <- popsfsSummary[popsfsSummary$Population=="African",]
pop_eur <- popsfsSummary[popsfsSummary$Population=="European",]
pop_easi <- popsfsSummary[popsfsSummary$Population=="East Asian",]

# SFS Plots by population ----
# African sfs
afr_plt = ggplot(data=sfs_afr, aes(x=AlleleFreq,y=mean,fill=ScalingFactor)) +
  theme_ponyo() +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), size=0.2, width=0.3, position=position_dodge(.9)) +
  ggtitle("Site Frequency Spectra of Human African Population") +
  labs(x="Minor Allele Frequency",y="Proportion of SNPs") + 
  scale_y_continuous(limits = c(0,0.65),expand = expansion(mult = c(0, .05))) +
  facet_nested(GenomeSize + DomCoefficent  ~ BurnInType) +
  theme(axis.text.x = element_text(angle = 50, vjust = 0.5))

# European sfs
eur_plt = ggplot(data=sfs_eur, aes(x=AlleleFreq,y=mean,fill=ScalingFactor)) +
  theme_ponyo() +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), size=0.2, width=0.3, position=position_dodge(.9)) +
  ggtitle("Site Frequency Spectra of Human European Population") +
  labs(x="Minor Allele Frequency",y="Proportion of SNPs") + 
  scale_y_continuous(limits = c(0,0.65),expand = expansion(mult = c(0, .05))) +
  facet_nested(GenomeSize + DomCoefficent  ~ BurnInType) +
  theme(axis.text.x = element_text(angle = 50, vjust = 0.5))

# East Asian sfs
easi_plt = ggplot(data=sfs_easi, aes(x=AlleleFreq,y=mean,fill=ScalingFactor)) +
  theme_ponyo() +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), size=0.2, width=0.3, position=position_dodge(.9)) +
  ggtitle("Site Frequency Spectra of Human East Asian Population") +
  labs(x="Minor Allele Frequency",y="Proportion of SNPs") + 
  scale_y_continuous(limits = c(0,0.65),expand = expansion(mult = c(0, .05))) +
  facet_nested(GenomeSize + DomCoefficent  ~ BurnInType) +
  theme(axis.text.x = element_text(angle = 50, vjust = 0.5))


# SFS tail plots by population ----
# African SFS Tail
tail_afr <- ggplot(data = sfsSummary[sfsSummary$Population=="African",], aes(x=BurnInType, y=TailSum_mean, fill=ScalingFactor))  + 
  theme_ponyo() +
  geom_errorbar(aes(ymin=TailSum_mean-TailSum_sd, ymax=TailSum_mean+TailSum_sd), width=.2, position=position_dodge(.9)) +
  ggtitle("SFS Tail Sums of Human African Population") +
  labs(x="Burn-In Method",y="Proportion of SNPs in SFS Tail\n(Sum of Allele Frequencies 5-50)") + 
  facet_grid2(DomCoefficent ~ GenomeSize, scales = "free_y", independent = "y")+
  scale_y_continuous(limits = c(0,0.6), expand = expansion(mult = c(0, .05)))

# European SFS Tail
tail_eur <- ggplot(data = sfsSummary[sfsSummary$Population=="European",], aes(x=BurnInType, y=TailSum_mean, fill=ScalingFactor))  + 
  theme_ponyo() +
  geom_errorbar(aes(ymin=TailSum_mean-TailSum_sd, ymax=TailSum_mean+TailSum_sd), width=.2, position=position_dodge(.9)) +
  ggtitle("SFS Tail Sums of Human European Population") +
  labs(x="Burn-In Method",y="Proportion of SNPs in SFS Tail\n(Sum of Allele Frequencies 5-50)") + 
  facet_grid2(DomCoefficent ~ GenomeSize, scales = "free_y", independent = "y")+
  scale_y_continuous(limits = c(0,0.6), expand = expansion(mult = c(0, .05)))

# East Asian SFS Tail
tail_easi <- ggplot(data = sfsSummary[sfsSummary$Population=="East Asian",], aes(x=BurnInType, y=TailSum_mean, fill=ScalingFactor))  + 
  theme_ponyo() +
  geom_errorbar(aes(ymin=TailSum_mean-TailSum_sd, ymax=TailSum_mean+TailSum_sd), width=.2, position=position_dodge(.9)) +
  ggtitle("SFS Tail Sums of Human East Asian Population") +
  labs(x="Burn-In Method",y="Proportion of SNPs in SFS Tail\n(Sum of Allele Frequencies 5-50)") + 
  facet_grid2(DomCoefficent ~ GenomeSize, scales = "free_y", independent = "y")+
  scale_y_continuous(limits = c(0,0.6), expand = expansion(mult = c(0, .05)))




# Heterozygosity by population ----

# African population heterozygosity
pi_afr <- ggplot(data = pop_afr, aes(x=BurnInType, y=ExpHet_mean, fill=ScalingFactor))  + 
  theme_ponyo() +
  geom_errorbar(aes(ymin=ExpHet_mean-ExpHet_sd, ymax=ExpHet_mean+ExpHet_sd), width=.2, position=position_dodge(.9)) +
  ggtitle("Heterozygosity of Human African Population") +
  labs(x="Coalescence Method",y="Expected Heterozygosity per 10kb") + 
  facet_grid(DomCoefficent ~ GenomeSize, scales = "free_y")+
  scale_y_continuous(limits =c(0,0.001), expand = expansion(mult = c(0, .05)),labels = label_number(suffix = "", scale = 10000)) +
  theme(axis.text.x = element_text(size = 12))

# European population heterozygosity
pi_eur <- ggplot(data = pop_eur, aes(x=BurnInType, y=ExpHet_mean, fill=ScalingFactor))  + 
  theme_ponyo() +
  geom_errorbar(aes(ymin=ExpHet_mean-ExpHet_sd, ymax=ExpHet_mean+ExpHet_sd), width=.2, position=position_dodge(.9)) +
  ggtitle("Heterozygosity of Human European Population") +
  labs(x="Coalescence Method",y="Expected Heterozygosity per 10kb") + 
  facet_grid(DomCoefficent ~ GenomeSize, scales = "free_y")+
  scale_y_continuous(limits =c(0,0.001), expand = expansion(mult = c(0, .05)),labels = label_number(suffix = "", scale = 10000)) +
  theme(axis.text.x = element_text(size = 12))

# East Asian population heterozygosity
pi_easi <- ggplot(data = pop_easi, aes(x=BurnInType, y=ExpHet_mean, fill=ScalingFactor))  + 
  theme_ponyo() +
  geom_errorbar(aes(ymin=ExpHet_mean-ExpHet_sd, ymax=ExpHet_mean+ExpHet_sd), width=.2, position=position_dodge(.9)) +
  ggtitle("Heterozygosity of Human East Asian Population") +
  labs(x="Coalescence Method",y="Expected Heterozygosity per 10kb") + 
  facet_grid(DomCoefficent ~ GenomeSize, scales = "free_y")+
  scale_y_continuous(limits =c(0,0.001), expand = expansion(mult = c(0, .05)),labels = label_number(suffix = "", scale = 10000)) +
  theme(axis.text.x = element_text(size = 12))

# Watterson's Theta by population ----

# African population Watterson's theta
theta_afr <- ggplot(data = pop_afr, aes(x=BurnInType, y=Theta_mean, fill=ScalingFactor))  + 
  theme_ponyo() +
  geom_errorbar(aes(ymin=Theta_mean-Theta_sd, ymax=Theta_mean+Theta_sd), width=.2, position=position_dodge(.9)) +
  ggtitle("Watterson's Theta of Human African Population") +
  labs(x="Coalescence Method",y="Watterson's Theta") + 
  facet_grid(DomCoefficent ~ GenomeSize, scales = "free_y") +
  scale_y_continuous(limits =c(0,0.0011), expand = expansion(mult = c(0, .05))) +
  theme(axis.text.x = element_text(size = 12))

# European population Watterson's theta
theta_eur <- ggplot(data = pop_eur, aes(x=BurnInType, y=Theta_mean, fill=ScalingFactor))  + 
  theme_ponyo() +
  geom_errorbar(aes(ymin=Theta_mean-Theta_sd, ymax=Theta_mean+Theta_sd), width=.2, position=position_dodge(.9)) +
  ggtitle("Watterson's Theta of Human European Population") +
  labs(x="Coalescence Method",y="Watterson's Theta") + 
  facet_grid(DomCoefficent ~ GenomeSize, scales = "free_y") +
  scale_y_continuous(limits =c(0,0.0011), expand = expansion(mult = c(0, .05))) +
  theme(axis.text.x = element_text(size = 12))

# East Asian population Watterson's theta
theta_easi <- ggplot(data = pop_easi, aes(x=BurnInType, y=Theta_mean, fill=ScalingFactor))  + 
  theme_ponyo() +
  geom_errorbar(aes(ymin=Theta_mean-Theta_sd, ymax=Theta_mean+Theta_sd), width=.2, position=position_dodge(.9)) +
  ggtitle("Watterson's Theta of Human East Asian Population") +
  labs(x="Coalescence Method",y="Watterson's Theta") + 
  facet_grid(DomCoefficent ~ GenomeSize, scales = "free_y") +
  scale_y_continuous(limits =c(0,0.0011), expand = expansion(mult = c(0, .05))) +
  theme(axis.text.x = element_text(size = 12))

# Save sfs plots ----
ggsave(afr_plt, filename = "C:/Users/tferrari/Desktop/SlimBenchmark/figures/sfs/african_sfs_WG.pdf", 
       device = cairo_pdf, width = 11, height = 8.5, units = "in")
ggsave(eur_plt, filename = "C:/Users/tferrari/Desktop/SlimBenchmark/figures/sfs/european_sfs_WG.pdf", 
       device = cairo_pdf, width = 11, height = 8.5, units = "in")
ggsave(easi_plt, filename = "C:/Users/tferrari/Desktop/SlimBenchmark/figures/sfs/eastasian_sfs.pdf", 
       device = cairo_pdf, width = 11, height = 8.5, units = "in")


# Save sfs tail plots ----
ggsave(tail_afr, filename = "C:/Users/tferrari/Desktop/SlimBenchmark/figures/sfs/african_sfsTail.pdf", 
       device = cairo_pdf, width = 11, height = 8.5, units = "in")
ggsave(tail_eur, filename = "C:/Users/tferrari/Desktop/SlimBenchmark/figures/sfs/european_sfsTail.pdf", 
       device = cairo_pdf, width = 11, height = 8.5, units = "in")
ggsave(tail_easi, filename = "C:/Users/tferrari/Desktop/SlimBenchmark/figures/sfs/eastasian_sfsTail.pdf", 
       device = cairo_pdf, width = 11, height = 8.5, units = "in")

# Save pop stat plots ----
ggsave(pi_afr, filename = "C:/Users/tferrari/Desktop/SlimBenchmark/figures/sfs/human_african_het_exonic.pdf", 
       device = cairo_pdf, width = 11, height = 8.5, units = "in")
ggsave(pi_eur, filename = "C:/Users/tferrari/Desktop/SlimBenchmark/figures/sfs/human_european_het_exonic.pdf", 
       device = cairo_pdf, width = 11, height = 8.5, units = "in")
ggsave(pi_easi, filename = "C:/Users/tferrari/Desktop/SlimBenchmark/figures/sfs/human_eastasian_het.pdf", 
       device = cairo_pdf, width = 11, height = 8.5, units = "in")

ggsave(theta_afr, filename = "C:/Users/tferrari/Desktop/SlimBenchmark/figures/sfs/human_african_theta_exonic.pdf", 
       device = cairo_pdf, width = 11, height = 8.5, units = "in")
ggsave(theta_eur, filename = "C:/Users/tferrari/Desktop/SlimBenchmark/figures/sfs/human_european_theta_exonic.pdf", 
       device = cairo_pdf, width = 11, height = 8.5, units = "in")
ggsave(theta_easi, filename = "C:/Users/tferrari/Desktop/SlimBenchmark/figures/sfs/human_eastasian_theta.pdf", 
       device = cairo_pdf, width = 11, height = 8.5, units = "in")
