# SETUP ----
# Load libraries
library(dplyr)
library(tidyverse)
library(car)
library(data.table)
options(scipen=0)

# Set working dir
setwd("C:/Users/tferrari/Desktop/SlimBenchmark/human/data")

# READ IN BURN COAL STATS ----

# Read-in coal stats from burn-in
burn_coal_data <- as.data.frame(fread("CoalStats.txt", sep="\t"))
burn_coal_data$BurnInType <- factor(as.character(burn_coal_data$BurnInType), levels = c("5","10","20"))
burn_coal_data$ScalingFactor <- factor(as.character(burn_coal_data$ScalingFactor), levels = c("1","5","10"))
burn_coal_data$GenomeSize <- factor(as.character(burn_coal_data$GenomeSize), levels = c("1e+05","1e+06","1e+07"))


# READ IN MAIN COAL STATS ----

main_coal_data <- as.data.frame(fread("main_coalescence_stats.txt", sep="\t"))

# Average 100kb coalstats
main_coal_data <- rbind.data.frame(main_coal_data[main_coal_data$GenomeSize==1e5,] %>% 
                                     mutate(RepNum=floor((main_coal_data[main_coal_data$GenomeSize==1e5,]$RepNum-1)/10)+1) %>%
                                     group_by(BurnInType,ScalingFactor,GenomeSize,DomCoefficent,BurnNum,RepNum) %>%
                                     summarise_at(as.character(c("#MissingCoalEvents","%MissingHeight","%MissingBranchLen")), list(mean=mean)) %>%
                                     rename_with(~ c("#MissingCoalEvents","%MissingHeight","%MissingBranchLen"), 7:9),
                                   main_coal_data[main_coal_data$GenomeSize!=1e5,])
colnames(main_coal_data)[6] <- "SimNum"

# Adjust
main_coal_data$BurnInType <- factor(as.character(main_coal_data$BurnInType), levels = c("5","10","20"))
main_coal_data$ScalingFactor <- factor(as.character(main_coal_data$ScalingFactor), levels = c("1","5","10"))
main_coal_data$GenomeSize <- factor(as.character(main_coal_data$GenomeSize), levels = c("1e+05","1e+06","1e+07"))


# READ IN SUMMARY STATS ----

# Read in simulated data
sfsData <- as.data.frame(read_delim("subsamp30_sfs_table.txt",delim="\t"))

# Combine 100kb genomes
sfsData <- rbind.data.frame(sfsData[sfsData$GenomeSize==1e5,] %>% 
                              mutate(SimNum=floor((sfsData[sfsData$GenomeSize==1e5,]$SimNum-1)/10)+1) %>%
                              group_by(Population,BurnInType,ScalingFactor,GenomeSize,DomCoefficent,BurnNum,SimNum) %>%
                              summarise_at(as.character(c("TotNumMuts",1:30)), list(sum=sum)) %>%
                              rename_with(~ c("TotNumMuts",1:30), 8:38),
                            sfsData[sfsData$GenomeSize!=1e5,])

# Expected heterozygosity function 
# takes in freqs vector of segregating sites, returns heterozygosity (2pq) UNAJUSTED BY GENSIZE
calc_het <- function(freqs, nSamp) {
  pq = (1:nSamp/(2*nSamp))*(1-1:nSamp/(2*nSamp))
  exp_het <- sum(2*pq*freqs)
  return(exp_het)
}

# Watterson's theta function
# takes in tot# of segregating sites, returns theta (K/a_n)) UNAJUSTED BY GENSIZE
calc_theta <-function(K,nSamp){
  a_n = sum(1/seq(1,(nSamp*2)-1,1))
  theta = K/a_n
  return(theta)
}

# Add popstats to dataframe
popsfsData <- sfsData %>% 
  rowwise() %>%
  mutate( ExpHet = calc_het(c_across(`1`:`30`),30)/GenomeSize, 
          Theta = calc_theta(TotNumMuts,30)/GenomeSize,
          Singletons = `1`/TotNumMuts) %>%
  select(Population:SimNum, ExpHet, Theta, Singletons)

# since 100kb are aggregated in sets of 10, divide by 10
popsfsData[popsfsData$GenomeSize==1e5,c("ExpHet","Theta")] = popsfsData[popsfsData$GenomeSize==1e5,c("ExpHet","Theta")]/10

# Filter and reformat data
popsfsData <- popsfsData[popsfsData$BurnInType %in% c(5,10,20),]
popsfsData$BurnInType <- factor(as.character(popsfsData$BurnInType), levels = c("5","10","20"))
popsfsData$ScalingFactor <- factor(as.character(popsfsData$ScalingFactor), levels = c("1","5","10"))
popsfsData$GenomeSize <- factor(as.character(popsfsData$GenomeSize), levels = c("1e+05","1e+06","1e+07"))

# READ IN LD STATS ----

# Read-in simulated data
ldData <- as.data.frame(read_delim("dist1000_10bins_ld_table.txt",delim="\t")) %>% 
  mutate(total_nSNPs = rowSums(select(., ends_with("_nSNPs"))),
         total_sumR2 = rowSums(select(., ends_with("_sumR2"))))

# combine 1e5 x 10 genomes
ldData <- rbind.data.frame(ldData[ldData$GenomeSize==1e5,] %>% 
                             mutate(SimNum=floor((ldData[ldData$GenomeSize==1e5,]$SimNum-1)/10)+1) %>%
                             group_by(Population,BurnInType,ScalingFactor,GenomeSize,DomCoefficent,BurnNum,SimNum) %>%
                             summarise_at(vars(total_nSNPs:total_sumR2), list(sum=sum)) %>%
                             rename_with(~ c("total_nSNPs","total_sumR2"), ends_with("_sum")),
                           ldData[ldData$GenomeSize!=1e5,c(colnames(ldData)[1:7],"total_nSNPs","total_sumR2")])

# Average R2 over # SNPs
ldData$avgR2 <- ldData$total_sumR2/ldData$total_nSNPs
ldData <- ldData %>% subset(select = -c(total_nSNPs,total_sumR2))

# Filter and reformat data
ldData <- ldData[ldData$BurnInType %in% c(5,10,20),]
ldData$BurnInType <- factor(as.character(ldData$BurnInType), levels = c("5","10","20"))
ldData$ScalingFactor <- factor(as.character(ldData$ScalingFactor), levels = c("1","5","10"))
#ldData$GenomeSize <- factor(as.character(ldData$GenomeSize), levels = c("1e+05","1e+06","1e+07"))
ldData$GenomeSize <- factor(as.character(ldData$GenomeSize), levels = c("1e+05","1e+06","1e+07"))


# COMPARE BURN COAL STATS ----

# Log scale axes plot
options(scipen=999)
logplot <- ggplot(burn_coal_data, aes(x=`#MissingCoalEvents`, y = `%MissingBranchLen`, color=BurnInType)) +
  geom_point(size = 4) +
  scale_x_log10() +
  scale_y_log10() +
  scale_colour_manual(values = c("5"="#E75B64" , 
                                 "10"= "#D8AF37", 
                                 "20"="#288B9A")) + 
  labs(x = "# Missing Coalescence Events", y="% Missing Branch Length") +
  theme_bw() + 
  theme(plot.title = element_text(size = 18, hjust = 0.5), 
        axis.text.x = element_text(size = 14, vjust = 1, hjust = 0.5), 
        axis.text.y = element_text(size = 14), 
        axis.title = element_text(size = 18),
        legend.title = element_text(size = 16), 
        legend.text = element_text(size = 16), 
        legend.position = "bottom")


# Linear regression fit
regfit<-lm(`%MissingBranchLen` ~ `#MissingCoalEvents`, data = burn_coal_data)

# Function to plot with linreg 
options(scipen=0)
ggplotRegression <- function(fit, color) {
  ggplot(fit$model, aes(x = .data[[names(fit$model)[2]]], y = .data[[names(fit$model)[1]]], color = color)) + 
    geom_point(size = 4) + 
    stat_smooth(method = 'lm', col = "black") +  
    theme_bw() +
    labs(title = bquote(R^2 == ~.(signif(summary(fit)$adj.r.squared, 4)) ~ "&" ~ "p" == ~.(signif(summary(fit)$coef[2,4], 4))))
}

# Plotting with updated theme and labels
regplot <- ggplotRegression(regfit, burn_coal_data$BurnInType) + 
  labs(x = "# Missing Coalescence Events", y="% Missing Branch Length") +
  scale_colour_manual(values = c("5"="#E75B64" , 
                                 "10"= "#D8AF37", 
                                 "20"="#288B9A")) +
  guides(color=guide_legend(title="Burn-in Length")) +
  theme_bw() +
  theme(plot.title = element_text(size = 18, hjust = 0.5), 
        axis.text.x = element_text(size = 14, vjust = 1, hjust = 0.5), 
        axis.text.y = element_text(size = 14), 
        axis.title = element_text(size = 18),
        legend.title = element_text(size = 16), 
        legend.text = element_text(size = 16), 
        legend.position = "bottom")

# MERGE MAIN COAL STATS AND SUMMARY STATS ----

# Merge popstat and LD data
popLDdata <- merge(ldData, popsfsData,
                   by = c("Population","BurnInType","ScalingFactor","GenomeSize","DomCoefficent","BurnNum","SimNum"),
                   all = TRUE)

# Merge coal stat and popstat data
main_coalpop <- merge(main_coal_data, popLDdata,
                      by = c("BurnInType","ScalingFactor","GenomeSize","DomCoefficent","BurnNum","SimNum"), 
                      all = TRUE)

# COAL VS EXPHET ----
options(scipen=999)

# Regression 
het_fit <- lm(ExpHet ~ `#MissingCoalEvents`*BurnInType + GenomeSize + ScalingFactor, data=main_coalpop[main_coalpop$Population=="african",])
summary(het_fit)

# Plot
coalVShet_afr <- ggplot(main_coalpop[main_coalpop$Population=="african",],
                        aes(x=`#MissingCoalEvents`, y = ExpHet, color=BurnInType)) +
  geom_point(size = 2) +
  scale_colour_manual(name = "Burn-In Length",
                      labels= c("5N", "10N","20N"),
                      values = c("5"="#E75B64" , 
                                 "10"= "#D8AF37", 
                                 "20"="#288B9A")) + 
  labs(x = "C-score", 
       y="Expected Heterozygosity") +
  theme_bw() + 
  theme(plot.title = element_text(size = 18, hjust = 0.5), 
        axis.text.x = element_text(size = 14, vjust = 1, hjust = 0.5), 
        axis.text.y = element_text(size = 14), 
        axis.title = element_text(size = 18),
        legend.title = element_text(size = 16), 
        legend.text = element_text(size = 16), 
        legend.position = "bottom")

# COAL VS THETA ----

# Regression 
theta_fit <- lm(Theta ~ `#MissingCoalEvents`*BurnInType + GenomeSize + ScalingFactor, data=main_coalpop[main_coalpop$Population=="african",])
summary(theta_fit)

# Plot
coalVStheta_afr <- ggplot(main_coalpop[main_coalpop$Population=="african",],
                        aes(x=`#MissingCoalEvents`, y = Theta, color=BurnInType)) +
  geom_point(size = 2) +
  scale_colour_manual(name = "Burn-In Length",
                      labels= c("5N", "10N","20N"),
                      values = c("5"="#E75B64" , 
                                 "10"= "#D8AF37", 
                                 "20"="#288B9A")) + 
  labs(x = "C-score", 
       y="Watterson's Theta") +
  theme_bw() + 
  theme(plot.title = element_text(size = 18, hjust = 0.5), 
        axis.text.x = element_text(size = 14, vjust = 1, hjust = 0.5), 
        axis.text.y = element_text(size = 14), 
        axis.title = element_text(size = 18),
        legend.title = element_text(size = 16), 
        legend.text = element_text(size = 16), 
        legend.position = "bottom")

# COAL VS LD ----

# Regression 
ld_fit <- lm(avgR2 ~ `#MissingCoalEvents`*BurnInType + GenomeSize + ScalingFactor, data=main_coalpop[main_coalpop$Population=="african",])
summary(ld_fit)

# Plot
coalVSld_afr <- ggplot(main_coalpop[main_coalpop$Population=="african",],
                          aes(x=`#MissingCoalEvents`, y = avgR2, color=BurnInType)) +
  geom_point(size = 2) +
  scale_colour_manual(name = "Burn-In Length",
                      labels= c("5N", "10N","20N"),
                      values = c("5"="#E75B64" , 
                                 "10"= "#D8AF37", 
                                 "20"="#288B9A")) + 
  labs(x = "C-score", 
       y="Average r²") + # Max 1Kb between SNPs
  theme_bw() + 
  theme(plot.title = element_text(size = 18, hjust = 0.5), 
        axis.text.x = element_text(size = 14, vjust = 1, hjust = 0.5), 
        axis.text.y = element_text(size = 14), 
        axis.title = element_text(size = 18),
        legend.title = element_text(size = 16), 
        legend.text = element_text(size = 16), 
        legend.position = "bottom")

# SAVE PLOTS ----

# Save compare coal stats plot
ggsave(logplot, filename = "C:/Users/tferrari/Desktop/SlimBenchmark/figures/Round2_GBE_Revisions/human_burnin_coaleventsVSbranchlen_log10axes.pdf",
       device = cairo_pdf, width = 11, height = 8.5, units = "in")
ggsave(regplot, filename = "C:/Users/tferrari/Desktop/SlimBenchmark/figures/Round2_GBE_Revisions/human_burnin_coaleventsVSbranchlen_linreg.pdf", 
       device = cairo_pdf, width = 11, height = 8.5, units = "in")

# Save coal vs summary stats
ggsave(coalVShet_afr, 
       filename = "C:/Users/tferrari/Desktop/SlimBenchmark/figures/GBE_Revisions_2/1_raw_figures/coalescence_events_vs/human_coalVShet.pdf", 
       device = cairo_pdf, width = 5, height = 4.25, units = "in")
ggsave(coalVStheta_afr, 
       filename = "C:/Users/tferrari/Desktop/SlimBenchmark/figures/GBE_Revisions_2/1_raw_figures/coalescence_events_vs/human_coalVStheta.pdf", 
       device = cairo_pdf, width = 5, height = 4.25, units = "in")
ggsave(coalVSld_afr, 
       filename = "C:/Users/tferrari/Desktop/SlimBenchmark/figures/GBE_Revisions_2/1_raw_figures/coalescence_events_vs/human_coalVSld.pdf", 
       device = cairo_pdf, width = 5, height = 4.25, units = "in")
