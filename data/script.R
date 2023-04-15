library("phyloseq")
library("ggplot2")
library("microbiome")
dyn.load('/usr/local/gfortran/lib/libgfortran.5.dylib')
library("DESeq2")
library("qiime2R")
library("data.table")
library("ggpubr")
library("vegan")
library("ape")
library("plyr")
library("RColorBrewer")
library("reshape2")
library("scales")
library("dplyr")
library("microbiomeutilities")
library("picante")




setwd("C:/Users/USER/Desktop/data/")


physeq<-qza_to_phyloseq(
  features="C:/Users/USER/Desktop/data/table-dada2.qza",
  tree="C:/Users/USER/Desktop/data/rooted-tree.qza",
  taxonomy="C:/Users/USER/Desktop/data/taxonomy.qza",
  metadata = "C:/Users/USER/Desktop/data/stats-dada2.txt"
)
ps <- physeq

###############################################################
## Alpha Diversity

GI <- microbiome::alpha(ps, index = "all")
sample_data(ps)$observed <- GI$observed
sample_data(ps)$Chao1 <- GI$chao1
sample_data(ps)$simpson <- GI$dominance_simpson
sample_data(ps)$Shannon <- GI$diversity_shannon

alpha_data <- sample_data(ps)
alpha_data <- alpha_data

theme_update(panel.background = element_rect(fill = "grey90"))
pdf("Alpha_diversity_all.pdf", width = 12)
##Observed_Diversity
ggplot(sample_data(ps), aes(x = Status, y = observed)) +
  ggtitle(paste("Observed_Diversity")) +
  geom_boxplot()+
  theme(plot.title = element_text(hjust = 0.5))+
  geom_point(aes( colour = Status))

##Cha01_Diversity 
ggplot(sample_data(ps), aes(x = Status, y = Chao1)) + 
  ggtitle(paste("Cha01_Diversity")) +
  geom_boxplot()+
  theme(plot.title = element_text(hjust = 0.5))+
  geom_point(aes( colour = Status))

##Shannon_Diversity
ggplot(sample_data(ps), aes(x = Status, y = Shannon)) + 
  ggtitle(paste("Shannon_Diversity")) +
  geom_boxplot()+
  theme(plot.title = element_text(hjust = 0.5))+
  geom_point(aes( colour = Status))

##simpson_Diversity
ggplot(sample_data(ps), aes(x = Status, y = simpson)) + 
  ggtitle(paste("simpson_Diversity")) +
  geom_boxplot()+
  theme(plot.title = element_text(hjust = 0.5))+
  geom_point(aes( colour = Status))

plot_richness(ps, color = "Status")+
  ggtitle(paste("Dervisty_by_Status")) +
  theme(plot.title = element_text(hjust = 0.5))

dev.off()


plot_richness(ps, color = "Status") +
  geom_point(size=3.5, alpha=0.7)



###############################################################
# Beta Diversity

carbom.ord <- ordinate(ps, "NMDS", "bray")
plot_ordination(ps, carbom.ord, type="taxa", color="Class", shape= "Kingdom", 
                title="OTUs")


plot_ordination(ps, carbom.ord, type="samples", color="Status", title="Status") + geom_point(size=3)

plot_ordination(ps, carbom.ord, type="taxa", color="Class", 
                title="OTUs", label="Class") + 
  facet_wrap(~Kingdom, 3)

###############################################################
## Taxa prevalence
ps <- subset_taxa(ps, !is.na(Species))
top20 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:20]
ps.top20 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
ps.top20 <- prune_taxa(top20, ps.top20)
plot_bar(ps.top20, fill="Phylum")



#####################################################################
## HeatMap
plot_heatmap(ps.top20, method = "NMDS", distance = "bray", 
             taxa.label = "Phylum", taxa.order = "Genus", 
             low="#FBFBF6F5", high="steelblue4", na.value="#FBFBF6F5")


#####################################################################
## Netwrok

total = median(sample_sums(ps))
standf = function(x, t=total) round(t * (x / sum(x)))
ps = transform_sample_counts(ps, standf)
ps_abund <- filter_taxa(ps, function(x) sum(x > total*0.10) > 0, TRUE)

plot_net(ps.top20, distance = "(A+B-2*J)/(A+B)", type = "taxa", 
         maxdist = 0.8, color="Kingdom", point_label="Phylum") 


#####################################################################
## LDA
library("microbiomeMarker")

mm_lefse <- run_lefse(
  ps,
  wilcoxon_cutoff = 0.5,
  group = "Status",
  kw_cutoff = 0.5,
  multigrp_strat = TRUE,
  lda_cutoff = 4
)

plot_cladogram(mm_lefse, color = c("darkgreen", "red"))

otu_table(ps)






