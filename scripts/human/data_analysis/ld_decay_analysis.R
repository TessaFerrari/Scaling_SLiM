# LOAD LIBRARIES AND READ IN SIMULATED DATA ----

# Load libraries
library(tidyverse)
library(extrafont)
library(scales)
library(ggpubr)
library(ggh4x)
library(ggbreak) 
library(assertr)
loadfonts(device = "win", quiet = TRUE)

# Set working directory
setwd("C:/Users/tferrari/Desktop/SlimBenchmark/human/data")

# Read-in simulated data
ldData <- as.data.frame(read_delim("dist1000_10bins_ld_table.txt",delim="\t"))


# READ IN EMPIRICAL DATA ----

multiIdx <- expand.grid(c("african"),c(5,10,20,"Coal","Recap"),c("Empirical"),c(1e5,1e6,1e7),c(0.0,0.5),c(1),c(1))
YRIdata <- as.data.frame(read_delim("empiricalData/dist1000_10bins_humanEmpirical_YRI.txt",delim="\t"))[,2:41] # Exonic data (totLength=26824075)
YRIdata <- rbind(YRIdata, YRIdata[rep(1, 29), ])
YRIdata <- cbind(multiIdx,YRIdata)  
colnames(YRIdata) <- c("Population", "BurnInType", "ScalingFactor", "GenomeSize", "DomCoefficent", "BurnNum", "SimNum", colnames(YRIdata)[8:47])


# COMBINE EVERY 10 100kb GENOMES INTO 1 1Mb GENOME ----

ldData <- rbind.data.frame(ldData[ldData$GenomeSize==1e5,] %>% 
                              mutate(SimNum=floor((ldData[ldData$GenomeSize==1e5,]$SimNum-1)/10)+1) %>%
                              group_by(Population,BurnInType,ScalingFactor,GenomeSize,DomCoefficent,BurnNum,SimNum) %>%
                              summarise_at(vars(b1_nSNPs:b10_Dp), list(sum=sum)) %>%
                              rename_with(~ col_concat(expand.grid(c("b"),c("_nSNPs","_sumR2","_sumD","_Dp"),as.character(1:10))[,c(1,3,2)]), ends_with("_sum")),
                            ldData[ldData$GenomeSize!=1e5,])


# COMBINE SIMULATED AND EMPIRICAL SFS ----

# Update factor levels
ldData$ScalingFactor <- as.character(ldData$ScalingFactor)

# Combine
ldData <- rbind(ldData,YRIdata)


# PIVOT BINS AND AVERAGE STATS OVER nSNPs ----

# Pivot
ldData_pivot <- ldData %>% 
  pivot_longer(cols = ends_with(c("_nSNPs","_sumR2","_sumD","_Dp")), 
               names_to = c("Bin", ".value"), 
               names_sep="_" ) 

# Average
ldData_pivot[,10:12] <- ldData_pivot[,10:12]/ldData_pivot$nSNPs


# SUMMARIZE STATS WITH MEAN AND SD----

ldStatSummary <- ldData_pivot %>%
  group_by(Population,BurnInType,ScalingFactor,GenomeSize,DomCoefficent,Bin) %>%
  summarise_at(c("sumR2","sumD","Dp"), 
               list(mean = mean, sd = sd))


# REFORMAT LD DATA ----

# Rename parameters and turn into factors
ldStatSummary$ScalingFactor <- factor(ldStatSummary$ScalingFactor, levels = c("Empirical","1","5","10"))
ldStatSummary$Population <- case_when(ldStatSummary$Population=="african" ~ "African",
                                   ldStatSummary$Population=="european" ~ "European",
                                   ldStatSummary$Population=="eastasian" ~ "East Asian")
ldStatSummary$BurnInType <- case_when(ldStatSummary$BurnInType=="5" ~ "5N",
                                   ldStatSummary$BurnInType=="10" ~ "10N",
                                   ldStatSummary$BurnInType=="20" ~ "20N",
                                   ldStatSummary$BurnInType=="Coal" ~ "Coal",
                                   ldStatSummary$BurnInType=="Recap" ~ "Recap")
ldStatSummary$BurnInType <- factor(ldStatSummary$BurnInType, levels=c('5N','10N','20N','Coal','Recap'))
ldStatSummary$GenomeSize <- case_when(ldStatSummary$GenomeSize==1e+05 ~ "100kb x 10",
                                   ldStatSummary$GenomeSize==1e+06 ~ "1Mb",
                                   ldStatSummary$GenomeSize==1e+07 ~ "10Mb")
ldStatSummary$GenomeSize <- factor(ldStatSummary$GenomeSize, levels=c('100kb x 10','1Mb','10Mb'))
ldStatSummary$DomCoefficent <- case_when(ldStatSummary$DomCoefficent==0.0 ~ "Recessive",
                                      ldStatSummary$DomCoefficent==0.5 ~ "Additive")
ldStatSummary$DomCoefficent <- factor(ldStatSummary$DomCoefficent, levels=c('Recessive','Additive'))

ldStatSummary$Bin <- case_when(ldStatSummary$Bin=="b1" ~ "0.0-0.1",
                               ldStatSummary$Bin=="b2" ~ "0.1-0.2",
                               ldStatSummary$Bin=="b3" ~ "0.2-0.3",
                               ldStatSummary$Bin=="b4" ~ "0.3-0.4",
                               ldStatSummary$Bin=="b5" ~ "0.4-0.5",
                               ldStatSummary$Bin=="b6" ~ "0.5-0.6",
                               ldStatSummary$Bin=="b7" ~ "0.6-0.7",
                               ldStatSummary$Bin=="b8" ~ "0.7-0.8",
                               ldStatSummary$Bin=="b9" ~ "0.8-0.9",
                               ldStatSummary$Bin=="b10" ~ "0.9-1.0")


# GGPLOT THEME FORMATTING FOR GRAPHS ----

# Reference colors
ponyoPalette= c("#4D413F","#5A7080","#288B9A","#E75B64","#DE7862","#D8AF37","#E8C4A2","#F8E7D3")

# Custom theme function
theme_ponyo <- function(){ 
  graphPalette= c("Empirical"="#5A7080","1"="#288B9A","5"="#E75B64","10"="#D8AF37")
  font = "Myriad Pro"   #assign font family up front
  dodge = position_dodge(width = 0.9)
  list(scale_fill_manual(name = "Scaling\nFactor", values = graphPalette),
       #scale_y_continuous(expand = expansion(mult = c(0, .05)),labels = function(x) format(x, scientific = TRUE)),
       geom_bar(stat = "identity", width=0.9, linewidth = 0.2, position = dodge, show.legend = TRUE, colour="black"),
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


# SEPARATE DATA BY POPULATION ----

ld_afr <- ldStatSummary[ldStatSummary$Population=="African",]
ld_eur <- ldStatSummary[ldStatSummary$Population=="European",]
ld_easi <- ldStatSummary[ldStatSummary$Population=="East Asian",]


# BINNED R2 PLOTS ----

# African R2
afr_r2_plt = ggplot(data=ld_afr, aes(x=Bin,y=sumR2_mean,fill=ScalingFactor)) +
  theme_ponyo() +
  geom_errorbar(aes(ymin=sumR2_mean-sumR2_sd, ymax=sumR2_mean+sumR2_sd), linewidth=0.2, width=0.3, position=position_dodge(.9)) +
  ggtitle("Coalescence Method") +
  labs(x="Distance between SNPs (Kb)",y="Average r²") + 
  scale_y_continuous(limits = c(0,0.6),expand = expansion(mult = c(0, .05)),
                     sec.axis = sec_axis(~ .,name = "Simulated Segment Length")) +
  facet_nested(GenomeSize + DomCoefficent  ~ BurnInType) +
  theme(axis.text.x = element_text(angle = 50, vjust = 0.5),
        axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        plot.title = element_text(size = 16, margin=margin(0,0,5,0)),
        legend.title = element_text(size = 10), legend.text = element_text(size = 8))

# Simple African R2
simp_ld_afr <- ld_afr[ld_afr$BurnInType=="Recap" & ld_afr$DomCoefficent=="Additive",]
simp_afr_r2_plt <- ggplot(data = simp_ld_afr, aes(x=Bin, y=sumR2_mean, fill=ScalingFactor))  + 
  theme_ponyo() +    #replace elements we want to change
  geom_errorbar(aes(ymin=sumR2_mean-sumR2_sd, ymax=sumR2_mean+sumR2_sd), size=0.2, width=0.3, position=position_dodge(.9)) +
  ggtitle("Simulated Segment Length") +
  labs(x="Distance between SNPs (Kb)",y="Average r²") + 
  facet_wrap(. ~ GenomeSize) +
  scale_y_continuous(limits = c(0,0.55),expand = expansion(mult = c(0, .05))) +
  theme(axis.text.x = element_text(angle = 50, vjust = 0.5),
        plot.title = element_text(size = 16))

# BINNED D PLOTS ----

# African D
afr_d_plt = ggplot(data=ld_afr, aes(x=Bin,y=sumD_mean,fill=ScalingFactor)) +
  theme_ponyo() +
  geom_errorbar(aes(ymin=sumD_mean-sumD_sd, ymax=sumD_mean+sumD_sd), linewidth=0.2, width=0.3, position=position_dodge(.9)) +
  ggtitle("Linkage Disequalibrium (D) of Human African Population") +
  labs(x="Distance between SNPs (Kb)",y="Average D") + 
  scale_y_continuous(limits = c(-0.04,0.09),expand = expansion(mult = c(.05, .05))) +
  facet_nested(GenomeSize + DomCoefficent  ~ BurnInType) +
  theme(axis.text.x = element_text(angle = 50, vjust = 0.5))


# BINNED D' PLOTS ----

# African D'
afr_dp_plt = ggplot(data=ld_afr, aes(x=Bin,y=Dp_mean,fill=ScalingFactor)) +
  theme_ponyo() +
  geom_errorbar(aes(ymin=Dp_mean-Dp_sd, ymax=Dp_mean+Dp_sd), linewidth=0.2, width=0.3, position=position_dodge(.9)) +
  ggtitle("Linkage Disequalibrium (D') of Human African Population") +
  labs(x="Distance between SNPs (Kb)",y="Average D'") + 
  scale_y_continuous(limits = c(-0.4,0.48),expand = expansion(mult = c(.05, .05))) +
  facet_nested(GenomeSize + DomCoefficent  ~ BurnInType) +
  theme(axis.text.x = element_text(angle = 50, vjust = 0.5))


# SAVE PLOTS TO PDF ----

# Save LD plots
ggsave(afr_r2_plt, filename = "C:/Users/tferrari/Desktop/SlimBenchmark/figures/GBE_Revisions_2/1_raw_figures/refurbished_figs/human_african_LD_r2_wEmpirical.pdf", 
       device = cairo_pdf, width = 11, height = 8.5, units = "in")
ggsave(simp_afr_r2_plt, filename = "C:/Users/tferrari/Desktop/SlimBenchmark/figures/GBE_Revisions_2/1_raw_figures/unfaceted_main_figs/human_african_LD_r2_wEmpirical_simple.pdf", 
       device = cairo_pdf, width = 11, height = 4.25, units = "in")
ggsave(afr_d_plt, filename = "C:/Users/tferrari/Desktop/SlimBenchmark/figures/LD/human_african_LD_D_wEmpirical.pdf", 
       device = cairo_pdf, width = 11, height = 8.5, units = "in")
ggsave(afr_dp_plt, filename = "C:/Users/tferrari/Desktop/SlimBenchmark/figures/LD/human_african_LD_Dprime_wEmpirical.pdf", 
       device = cairo_pdf, width = 11, height = 8.5, units = "in")

