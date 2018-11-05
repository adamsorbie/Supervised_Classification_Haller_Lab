library(dplyr)
library(ggplot2)

alpha_div_cc <- read.table("cc/alpha-diversity.tab", sep="\t", header = T, row.names = 1)
alpha_div_tissue <- read.table("tissue/alpha-diversity.tab", sep = "\t", header = T, row.names = 1)
mapping_cc <- read.table("cc/mapping_file_wo_controls.tab", sep = "\t", header = T, row.names = 1, comment.char = "")
mapping_tissue <- read.table("tissue/mapping_file_wo_controls.tab", sep = "\t", header = T, row.names = 1, comment.char = "")
## combine files 
combined_cc <- merge(alpha_div_cc, mapping_cc, by = 0)
row.names(combined_cc) <- combined_cc$Row.names
combined_cc <- subset(combined_cc, select = -c(Row.names) )
combined_cc_T <- combined_cc[combined_cc$Phenotype == "T",]

combined_tissue <- merge(alpha_div_tissue, mapping_tissue, by = 0)
row.names(combined_tissue) <- combined_tissue$Row.names
combined_tissue <- subset(combined_tissue, select = -c(Row.names))
combined_tissue_T <- combined_tissue[combined_tissue$Phenotype == "T",]


shannon_grouped_mean_cc <- aggregate(combined_cc$Shannon.effective, by=list(Age=combined_cc$Age, 
                                                                               Genotype=combined_cc$Genotype), mean)

shannon_grouped_mean_cc$site <- c("caecal")

shannon_grouped_mean_tissue <- aggregate(combined_tissue$Shannon.effective, by=list(Age=combined_tissue$Age,
                                                                                      Genotype=combined_tissue$Genotype), mean)

shannon_grouped_mean_tissue$site <- c("tissue")

merged_shannon <- rbind(shannon_grouped_mean_cc, shannon_grouped_mean_tissue)
genotypes_avi <- c("tg/cre;-/-", "tg/wt;-/-")
merged_shannon$Mouse_Line <- ifelse(merged_shannon$Genotype %in% genotypes_avi, "AVI", "AV")
merged_shannon_no_tgwt <- merged_shannon[merged_shannon$Genotype != "tg/wt",]


#shannon_grouped_mean <- shannon_grouped_mean[order(shannon_grouped_mean$Age),]
#shannon_grouped_mean["Phenotype"] <- c("NT", "NT", "NT", "NT", "T", "NT", "NT", "T", "NT" )

shannon_plot <- ggplot(data = merged_shannon_no_tgwt, aes(x = Age , y= x, group= Genotype, color=Genotype)) +
  geom_point(size=4) +
  geom_line() +
  facet_grid(merged_shannon_no_tgwt$site~merged_shannon_no_tgwt$Mouse_Line) +
  theme_bw() +
  xlab("Age (weeks)") +
  ylab("Shannon Diversity") +
  ggtitle("Shannon diversity across time") + theme(plot.title = element_text(hjust = 0.5))


shannon_plot + scale_x_discrete(limits=c("5wk","12wk","20wk")) 