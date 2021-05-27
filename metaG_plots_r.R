# R script to generate plots and statistics from metagenomic recruitment RPKM values

# This is meant to be run in R studio interactively and has different functions 
# commented out depending on which datasets one is analyzing. For example,
# the first two read.csv commands toggle between the global (Tara, Biogeotraces, Malaspina)
# and coastal (Dead Zone, SF Bay, So Cal Bight) RPKM datasets.
# Similar functionality for saving plots and alternating between different subsets of data
# occur throughout.


# Clear existing settings and data
rm(list=ls()) 

# Required libraries
require(ComplexHeatmap)
require(reshape2)
require(ggplot2)
require(maps)
#require(tidyverse)
#require(RColorBrewer)

##############################################################
# Read in the data 
# For global data:
df <- read.csv("~/Desktop/OM252_all_metaG_data.csv", header=T)  

# For coastal data:
df <- read.csv("~/Desktop/OM252_all_coastal_metaG_data_r.csv", header=T)  

# Transform matrix into ggplot-friendly format
datamelt=melt(df, id.vars=c("db", "sampling_date", "sampling_time", "sample", "thrash_db_acc", "depth", "depth_zone", "lat", "lon", 
                          "true_lat", "region", "sal", "sal_group", "temp", "temp_group"), 
          measure.vars=c(	"UBA12265", "UBA9605", "UBA9601", "UBA11194", "AG.896.J04", "UBA8357", "TOBG_MED.626", "AG.898.O07", 
                            "AG.915.K04",	"AG.359.J14", "AG.900.B21", "AG.918.E15", "TOBG_MED.814", "HIMB30", "AG.905.C17", "LSUCC0096", 
                            "TOBG_MED.759", "TOBG_RS.469", "TOBG_SAT.133", "TOBG_NAT.109", "GCA_002480175", "TOBG_NP.1444", 
                            "TOBG_SP.353", "UBA11144", "TOBG_NP.1472"))
colnames(datamelt)[16] <- "genome"
colnames(datamelt)[17] <- "RPKM"

# Set the arrangement order of the genomes according to the phylogenetic tree
datamelt$genome=factor(datamelt$genome, levels=c("AG.359.J14", "AG.900.B21", "AG.918.E15", "TOBG_MED.814", "AG.898.O07", "AG.915.K04", 
                                          "HIMB30", "AG.905.C17", "TOBG_MED.759", "LSUCC0096", "UBA12265", "UBA9605", "UBA9601", 
                                          "UBA11194", "TOBG_MED.626", "AG.896.J04", "UBA8357", "TOBG_RS.469", "TOBG_SAT.133", 
                                          "TOBG_NAT.109", "GCA_002480175", "TOBG_NP.1444", "TOBG_SP.353", "UBA11144", "TOBG_NP.1472"))

##############################################################
# Plot RPKM values according to site location on a world map
# This was generated through help from Murat Eren and his script 
# here: https://github.com/merenlab/world-map-r/

# Constrain limits of the map based on plotted values
min_lat <- min(datamelt$lat) + -5
max_lat <- max(datamelt$lat) + 5
min_lon <- min(datamelt$lon) + -5
max_lon <- max(datamelt$lon) + 5

# Generate the world plot
world <- map_data("world")
gg1 <- ggplot() +
  geom_polygon(data=world, aes(x=long, y=lat, group=group)) +
  coord_fixed(1.2, xlim=c(min_lon, max_lon), ylim=c(min_lat, max_lat)) +
  #theme_minimal()
  theme_classic()

# Add the RPKM data
gg2 <- gg1 + geom_jitter(data=datamelt, aes_string(x="lon", y="lat", size="RPKM", col="genome"), alpha=0.6, width=1.0, height=1.0) +
  scale_size_continuous(limits=c(0.01,3.0), breaks=c(0.1,0.5,1.0,2.0,3.0)) +
  labs(x="Longitude", y="Latitude", size="log10 RPKM")

# Plot
gg2


# To see the genome data on separate maps by genome
gg3 <- gg1 + geom_point(data=datamelt, aes_string(x="lon", y="lat", size="RPKM"), col="red", alpha=0.4) +
  #scale_size_continuous(limits=c(0.01,3.0), breaks=c(0.1,0.5,1.0,2.0,3.0)) +
  scale_radius(range=c(0.1,3)) +
  labs(x="Longitude", y="Latitude", size="log10 RPKM") +
  theme_classic()

# Repeat by genome
gg4 <- gg3 + facet_wrap(~ genome)

# Plot
gg4

ggsave("map_RPKM_by_genome.pdf", dpi = 300, width = 11, height = 8.5)

##############################################################
# Map abundances for specific regions
sub_data <- datamelt[datamelt$db=='SIERADZKI',]
#sub_data <- datamelt[datamelt$db=='THRASH',]

# Constrain limits of the map based on plotted values
# Adjust by location to achieve best visualization
min_lat <- min(sub_data$lat) + -0.5
max_lat <- max(sub_data$lat) + 0.5
min_lon <- min(sub_data$lon) + -0.5
max_lon <- max(sub_data$lon) + 0.5

# Generate the spatial plot
world <- map_data("world")
gg1 <- ggplot() +
  geom_polygon(data=world, aes(x=long, y=lat, group=group)) +
  coord_fixed(1.2, xlim=c(min_lon, max_lon), ylim=c(min_lat, max_lat)) +
  #theme_minimal()
  theme_classic()

# Add the RPKM data
gg2 <- gg1 + geom_jitter(data=sub_data, aes_string(x="lon", y="lat", size="RPKM", col="genome"), alpha=0.6, width=0.01, height=0.01) +
  scale_size_continuous(limits=c(0.01,3.0), breaks=c(0.1,0.5,1.0,2.0,3.0)) +
  labs(x="Longitude", y="Latitude", size="log10 RPKM")

# Plot
gg2


# To see the genome data on separate maps by genome
gg3 <- gg1 + geom_point(data=sub_data, aes_string(x="lon", y="lat", size="RPKM"), col="red", alpha=0.4) +
  #scale_size_continuous(limits=c(0.01,3.0), breaks=c(0.1,0.5,1.0,2.0,3.0)) +
  scale_radius(range=c(0.1,3)) +
  labs(x="Longitude", y="Latitude", size="log10 RPKM") +
  theme_classic() +
  theme(axis.text.x = element_text(angle=90, size=7))

# Repeat by genome
gg4 <- gg3 + facet_wrap(~ genome)

# Plot
gg4

ggsave("LA_map_RPKM_by_genome.pdf", dpi = 300, width = 11, height = 8.5)
#ggsave("DZ_map_RPKM_by_genome.pdf", dpi = 300, width = 11, height = 8.5)

##############################################################
# Plot stacked bar blot of all genome abundances by site
ggplot(datamelt) + 
  geom_bar(aes(x=sample, y=RPKM, fill=genome), stat="identity", position="stack") +
  labs(x="Samples", y="log10 RPKM") + 
  labs(fill="genome") +
  #scale_fill_manual(values = c("springgreen4","dodgerblue4","skyblue2")) +
  theme_classic() +
  theme(axis.text.x = element_text(angle=90, size=5))

# Plot side-by-side bar blot of all genome abundances by all sites
# or subset the data to examine specific sites

sub_data <- datamelt[datamelt$db=='SIERADZKI',]
#sub_data <- datamelt[datamelt$db=='THRASH',]

#ggplot(datamelt) + 
ggplot(sub_data) +
  geom_bar(aes(x=sample, y=RPKM, fill=genome), stat="identity", position=position_dodge()) +
  labs(x="Samples", y="log10 RPKM") + 
  labs(fill="genome") +
  theme_classic() +
  theme(axis.text.x=element_text(angle=90, size=8))


##############################################################
# Plot histograms and boxplots of all genomes separately by RPKM to explore relative abundances

ggplot(data=datamelt, aes(x=RPKM)) +
  geom_histogram(binwidth=0.1) +
  labs(x="log10 RPKM", y="Number") +
  facet_wrap(~ genome) + 
  theme_classic()

ggplot() +
  geom_boxplot(data=datamelt, aes(y=RPKM, x=genome), outlier.shape = 1) +
  labs(x="Genome", y="log10 RPKM") +
  theme_classic() +
  theme(axis.text.x=element_text(angle=90, size=8))

ggsave("genomes_by_rpkm_boxplot.pdf", dpi = 300, width = 11, height = 8.5)
#ggsave("genomes_by_rpkm_boxplot_coastal.pdf", dpi = 300, width = 11, height = 8.5)

# Subset for one of the coastal datasets
sub_data <- datamelt[datamelt$db=='SIERADZKI',]
#sub_data <- datamelt[datamelt$db=='THRASH',]

ggplot() +
  geom_boxplot(data=sub_data, aes(y=RPKM, x=genome), outlier.shape = 1) +
  labs(x="Genome", y="log10 RPKM") +
  theme_classic() +
  theme(axis.text.x=element_text(angle=90, size=8))


##############################################################
# Plot RPKM values of all genomes together according to site salinity, temperature, and depth
# and plot a linear trendline with 95% confidence intervals

# Temperature
ggplot() +
  geom_jitter(data=datamelt, aes_string(x="RPKM", y="temp", col="genome"), alpha=0.4, width=0.05, height=0.05) +
  geom_smooth(data=datamelt, aes(x=RPKM, y=temp), method='lm', color="grey") +
  labs(x="log10 RPKM", y="Temp (˚C)") +
  theme_classic()

# Salinity
ggplot() +
  geom_jitter(data=datamelt, aes_string(x="RPKM", y="sal", col="genome"), alpha=0.4, width=0.05, height=0.05) +
  geom_smooth(data=datamelt, aes(x=RPKM, y=sal), method='lm', color="grey") +
  labs(x="RPKM", y="Salinity") +
  theme_classic()

# Depth
ggplot() +
  geom_jitter(data=datamelt, aes_string(x="RPKM", y="depth", col="genome"), alpha=0.4, width=0.05, height=0.05) +
  scale_y_reverse(lim=c(6000,0)) +
  geom_smooth(data=datamelt, aes(x=RPKM, y=depth), method='lm', color="grey") +
  labs(x="RPKM", y="Depth (m)") +
  theme_classic()


##############################################################
# Separate plots by genome for RPKM values according to site salinity, temperature, and depth
# with a linear trendline with 95% confidence intervals

# Plot OM252 genome RPKMs vs temperature
ggplot() +
  geom_point(data=datamelt, aes(x=RPKM, y=temp), color="red", alpha=0.4, size=0.5) +
  labs(x="log10 RPKM", y="Temp (˚C)") +
  geom_smooth(data=datamelt, aes(x=RPKM, y=temp), method='lm', color="grey") +
  facet_wrap(~ genome, scales="free_x") + 
  theme_classic()

ggsave("Temp_v_RPKM.pdf", dpi = 300, width = 11, height = 8.5)

# Evaluate R2 and P values for the linear trends with different selected genomes:
# help from https://www.scribbr.com/statistics/linear-regression-in-r/
#sub_data_LSUCC0096 <- datamelt[datamelt$genome=='LSUCC0096',]
#plot(sal ~ RPKM, data=sub_data_LSUCC0096)
#relationship <- lm(sal ~ RPKM, data=sub_data_LSUCC0096)
#summary(relationship)
#par(mfrow=c(2,2))
#plot(relationship)

#sub_data_HIMB30 <- datamelt[datamelt$genome=='HIMB30',]
#plot(sal ~ RPKM, data=sub_data_HIMB30)
#relationship <- lm(sal ~ RPKM, data=sub_data_HIMB30)
#summary(relationship)
#par(mfrow=c(2,2))
#plot(relationship)
  
# Plot OM252 genome RPKMs vs salinity
ggplot() +
  geom_point(data=datamelt, aes(x=RPKM, y=sal), color="red", alpha=0.4, size=0.5) +
  labs(x="log10 RPKM", y="Salinity") +
  geom_smooth(data=datamelt, aes(x=RPKM, y=sal), method='lm', color="grey") +
  facet_wrap(~ genome, scales="free_x") + 
  theme_classic()

ggsave("Sal_v_RPKM.pdf", dpi = 300, width = 11, height = 8.5)
#ggsave("coastal_Sal_v_RPKM.pdf", dpi = 300, width = 11, height = 8.5)

# Evaluate R2 and P values for the linear trends with different selected genomes:
# help from https://www.scribbr.com/statistics/linear-regression-in-r/
#sub_data_LSUCC0096 <- datamelt[datamelt$genome=='LSUCC0096',]
#plot(sal ~ RPKM, data=sub_data_LSUCC0096)
#relationship <- lm(sal ~ RPKM, data=sub_data_LSUCC0096)
#summary(relationship)
#par(mfrow=c(2,2))
#plot(relationship)

#sub_data_HIMB30 <- datamelt[datamelt$genome=='HIMB30',]
#plot(sal ~ RPKM, data=sub_data_HIMB30)
#relationship <- lm(sal ~ RPKM, data=sub_data_HIMB30)
#summary(relationship)
#par(mfrow=c(2,2))
#plot(relationship)

# Plot OM252 genome RPKMs vs depth with non-linear trendline
ggplot() +
  geom_point(data=datamelt, aes(x=RPKM, y=depth), color="red", alpha=0.4, size=0.5) +
  scale_y_reverse(lim=c(6000,0)) +
  labs(x="log10 RPKM", y="Depth (m)") +
  geom_smooth(data=datamelt, aes(x=RPKM, y=depth), method='auto', color="grey") +
  facet_wrap(~ genome, scales="free_x") + 
  theme_classic()

ggsave("Depth_v_RPKM_nonlinear.pdf", dpi = 300, width = 11, height = 8.5)

# Plot OM252 genome RPKMs vs depth with linear trendline
ggplot() +
  geom_point(data=datamelt, aes(x=RPKM, y=depth), color="red", alpha=0.4, size=0.5) +
  scale_y_reverse(lim=c(6000,0)) +
  labs(x="log10 RPKM", y="Depth (m)") +
  geom_smooth(data=datamelt, aes(x=RPKM, y=depth), method='lm', color="grey") +
  facet_wrap(~ genome, scales="free_x") + 
  theme_classic()

ggsave("Depth_v_RPKM_linear.pdf", dpi = 300, width = 11, height = 8.5)

##############################################################
# Plot RPKM values according to latitude 

# Plot all together
ggplot() +
  geom_jitter(data=datamelt, aes_string(x="lat", y="RPKM", col="genome"), alpha=0.4, width=0.05, height=0.05) +
  labs(x="Latitude", y="log10 RPKM") +
  geom_smooth(data=datamelt, aes(x=lat, y=RPKM), method="auto", color="grey") +
  theme_classic()

# Plot each genome separately 
ggplot() +
  geom_point(data=datamelt, aes(x=lat, y=RPKM), color="red", alpha=0.4, size=0.5) +
  labs(x="Latitude", y="log10 RPKM") +
  geom_smooth(data=datamelt, aes(x=lat, y=RPKM), method = "auto", color="grey") +
  facet_wrap(~ genome, scales="free_y") + 
  theme_classic()

ggsave("Lat_v_RPKM.pdf", dpi = 300, width = 11, height = 8.5)

# Check sampling site location distribution...
ggplot(data=df, aes(x=lat, fill=db)) +
  geom_histogram(binwidth=1) +
  labs(x="Latitude", y="number of sampling sites", fill="dataset") +
  theme_classic()

ggsave("numsites_v_lat.pdf", dpi = 300, width = 11, height = 8.5)

# ...and density
ggplot(data=df, aes(x=lat)) +
  geom_histogram(aes(y=..density..), position="identity", alpha=1, binwidth=1) +
  geom_density(alpha=0.2, fill="#FF6666") +
  labs(x="Latitude", y="density of sampling sites") +
  theme_classic()

ggsave("densites_v_lat.pdf", dpi = 300, width = 11, height = 8.5)

##############################################################
# Plot latitude by genome via "heatmap" (geom_tile)

ggplot(datamelt, aes(x=lat, y=variable, fill=value)) +
  geom_tile() +
  #theme_minimal() +
  scale_fill_gradient(low="grey", high="red", na.value="white") +
  theme(axis.text.x = element_text(angle=90, size=5))

