library(reshape2)
library(ggplot2)
library(ggpubr)


alpha_div <- read.table("alpha-diversity.tab", sep="\t", header = T, row.names = 1, check.names = F, comment.char = "")
mapping <- read.table("mapping_file.tab", sep="\t", header = T, row.names = 1, check.names = F, comment.char = "")
combined <- merge(alpha_div, mapping, by = 0)
row.names(combined) <- combined$Row.names
combined <- subset(combined, select = -c(Row.names) )

# this code would also be much cleaner if you write plotting as a function 

# function to select column for plotting and melt df 

select_var <- function(df, column){
  cols <- c("Richness", "Shannon.effective", column)
  df1 <- df[cols]
  df1_m <- melt(df1)
  return(df1_m)
}

combined_12_20_only <- combined[combined$Age != "5wk", ] # plot 12 and 20wk only for R NR comparison
combined_5_only <- combined[combined$Age == "5wk", ] # plot 5wk only for 5wk genotype analysis

alpha_gen_a <- select_var(combined, "Genotype_Age")
alpha_gen_richness <- alpha_gen_a[alpha_gen_a$variable == "Richness", ]
alpha_gen_shannon.e <- alpha_gen_a[alpha_gen_a$variable == "Shannon.effective", ]

alpha_RS <- select_var(combined_12_20_only, "Responder")
alpha_RS_richness <- alpha_RS[alpha_RS$variable == "Richness", ]
alpha_RS_shannon.e <- alpha_RS[alpha_RS$variable == "Shannon.effective", ]

alpha_tissue_RS <- select_var(combined_12_20_only, "Tissue_Responder")
alpha_tissue_RS_richness <- alpha_tissue_RS[alpha_tissue_RS$variable == "Richness", ]
alpha_tissue_RS_shannon.e <- alpha_tissue_RS[alpha_tissue_RS$variable == "Shannon.effective", ]

alpha_tissue <- select_var(combined_12_20_only, "Tissue")
alpha_tissue_richness <- alpha_tissue[alpha_tissue$variable == "Richness", ]
alpha_tissue_shannon.e <- alpha_tissue[alpha_tissue$variable == "Shannon.effective", ]

alpha_gen_5 <- select_var(combined_5_only, "Genotype")
alpha_gen__5_richness <- alpha_gen_5[alpha_gen_5$variable == "Richness", ]
alpha_gen__5_shannon.e <- alpha_gen_5[alpha_gen_5$variable == "Shannon.effective", ]

######### plotting function 

plot_alpha <- function(df_richness,df_shannon.e, col, comparisons_list, stat_method){
                       compare_means(value ~ df_richness[[col]], data = df_richness, method = stat_method)
                       plot_richness <- ggplot(df_richness, aes(x=df_richness[[col]], y=df_richness["value"], 
                                               fill = df_richness[[col]] )) + labs(x=col, y="alpha diversity") + 
                       stat_boxplot(geom = "errorbar", width=0.2)
                       final_plot_richness <- plot_richness + stat_compare_means(comparisons = comparisons_list, label = "p.signif" ) + 
                       theme_bw() + ggtitle("Richness") + theme(axis.text.x = element_text(size = 8), axis.line = element_line(colour = "black"), 
                       panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),
                       panel.border = element_blank(), plot.title(element_text(hjust = 0.5))) 
                       
                       compare_means(value ~ df_shannon.e[[col]], data = df_shannon.e, method = stat_method)
                       plot_shannon.e <- ggplot(df_shannon.e, aes(x=df_shannon.e[[col]], y=df_shannon.e["value"], 
                                                                 fill = df_shannon.e[[col]] )) + labs(x=col, y="alpha diversity") + 
                       stat_boxplot(geom = "errorbar", width=0.2)
                       final_plot_shannon.e <- plot_shannon.e + stat_compare_means(comparisons = comparisons_list, label = "p.signif" ) + 
                       theme_bw() + ggtitle("Shannon.effective") + theme(axis.text.x = element_text(size = 8), axis.line = element_line(colour = "black"), 
                       panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),
                       panel.border = element_blank(), plot.title(element_text(hjust = 0.5)))
                       
                       return(final_plot_richness, final_plot_shannon.e)
}

####################################################### TISSUE_RESPONDER ####################################################################
compare_means(value ~ Tissue_Responder, data = alpha_tissue_RS_richness, method = "kruskal.test")
compare_means(value ~ Tissue_Responder, data = alpha_tissue_RS_shannon.e, method = "kruskal.test")

my_comparisons_RS <- list( c("Healthy_WT", "Healthy_NR"), c("Healthy_WT", "Tumor_R"), c("Healthy_WT", "Tumor_adjacent_R"), 
                        c("Tumor_R", "Tumor_adjacent_R"), c("Tumor_R", "Healthy_NR"), c("Tumor_adjacent_R", "Healthy_NR") )


p_richness_RS <- ggplot(alpha_tissue_RS_richness , aes(x=Tissue_Responder, y=value, fill= Tissue_Responder )) + geom_boxplot(color="black") +  
                     labs(x="Sample", y="alpha diversity") + 
                     stat_boxplot(geom ='errorbar', width=0.2) 

p_richness_RS + stat_compare_means(comparisons = my_comparisons_RS, label = "p.signif" ) + theme_bw() + 
                                   theme(axis.text.x = element_text(size = 8)) + ggtitle("Richness") + 
                                   theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(),
                                   panel.grid.minor = element_blank(), panel.border = element_blank(),
                                   panel.background = element_blank(), plot.title = element_text(hjust = 0.5))

p_shannon.e_RS <-  ggplot(alpha_tissue_RS_shannon.e , aes(x=Tissue_Responder, y=value, fill= Tissue_Responder)) + geom_boxplot(color="black") + 
                     labs(x="Sample", y="alpha diversity") + 
                     stat_boxplot(geom ='errorbar', width = 0.2)  
p_shannon.e_RS + stat_compare_means(comparisons = my_comparisons_RS, label ="p.signif") + theme_bw() + ggtitle("Shannon Effective") + 
                               theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(),
                               panel.grid.minor = element_blank(), panel.border = element_blank(),
                               panel.background = element_blank(), plot.title = element_text(hjust = 0.5)) + 
                               theme(axis.text.x = element_text(size = 8))

###############################################    GENOTYPE_AGE    ############################################################################

compare_means(value ~ Genotype_Age, data = alpha_gen_richness, method = "kruskal.test")
compare_means(value ~ Genotype_Age, data = alpha_gen_shannon.e, method = "kruskal.test")

my_comparisons_ga <- list( c("wt_5wk", "wt_12wk"), c("wt_12wk", "wt_20wk"), c("wt_5wk", "cre_5wk"), 
                           c("wt_12wk", "cre_12wk"), c("wt_20wk", "cre_20wk"), c("cre_5wk", "cre_12wk"),  
                           c("cre_5wk", "cre_20wk"), c("cre_12wk", "cre_20wk"))


p_richness_ga <- ggplot(alpha_gen_richness , aes(x=Genotype_Age, y=value, fill= Genotype_Age )) + geom_boxplot(color="black") +  
  labs(x="Sample", y="alpha diversity") + 
  stat_boxplot(geom ='errorbar', width=0.2) 

p_richness_ga + stat_compare_means(comparisons = my_comparisons_ga, label = "p.signif" ) + theme_bw() + 
  theme(axis.text.x = element_text(size = 8)) + ggtitle("Richness") + 
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), panel.border = element_blank(),
        panel.background = element_blank(), plot.title = element_text(hjust = 0.5))

p_shannon.e_ga <-  ggplot(alpha_gen_shannon.e , aes(x=Genotype_Age, y=value, fill= Genotype_Age)) + geom_boxplot(color="black") + 
  labs(x="Sample", y="alpha diversity") + 
  stat_boxplot(geom ='errorbar', width = 0.2)  
p_shannon.e_ga + stat_compare_means(comparisons = my_comparisons_ga, label ="p.signif") + theme_bw() + ggtitle("Shannon Effective") + 
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), panel.border = element_blank(),
        panel.background = element_blank(), plot.title = element_text(hjust = 0.5)) + 
  theme(axis.text.x = element_text(size = 8))


############################################################## TISSUE #########################################################################

compare_means(value ~ Tissue, data = alpha_tissue_richness, method = "kruskal.test")
compare_means(value ~ Tissue, data = alpha_tissue_shannon.e, method = "kruskal.test")

my_comparisons_tissue <- list( c("Healthy", "Tumor"), c("Healthy", "Tumor_adjacent"), c("Tumor", "Tumor_adjacent"))


p_richness_tissue <- ggplot(alpha_tissue_richness , aes(x=Tissue, y=value, fill= Tissue )) + geom_boxplot(color="black") +  
  labs(x="Tissue", y="alpha diversity") + 
  stat_boxplot(geom ='errorbar', width=0.2) 

p_richness_tissue + stat_compare_means(comparisons = my_comparisons_tissue, label = "p.signif" ) + theme_bw() + 
  theme(axis.text.x = element_text(size = 8)) + ggtitle("Richness") + 
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), panel.border = element_blank(),
        panel.background = element_blank(), plot.title = element_text(hjust = 0.5))

p_shannon.e_tissue <-  ggplot(alpha_tissue_shannon.e , aes(x=Tissue, y=value, fill= Tissue)) + geom_boxplot(color="black") + 
  labs(x="Tissue", y="alpha diversity") + 
  stat_boxplot(geom ='errorbar', width = 0.2)  
p_shannon.e_tissue + stat_compare_means(comparisons = my_comparisons_tissue, label ="p.signif") + theme_bw() + ggtitle("Shannon Effective") + 
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), panel.border = element_blank(),
        panel.background = element_blank(), plot.title = element_text(hjust = 0.5)) + 
  theme(axis.text.x = element_text(size = 8))

############################################################## RESPONDER #########################################################################

compare_means(value ~ Responder, data = alpha_RS_richness, method = "kruskal.test")
compare_means(value ~ Responder, data = alpha_RS_shannon.e, method = "kruskal.test")

my_comparisons_RS <- list( c("R", "NR"), c("WT", "R"), c("WT", "NR"))


p_richness_RS <- ggplot(alpha_RS_richness , aes(x=Responder, y=value, fill= Responder )) + geom_boxplot(color="black") +  
  labs(x="Responder", y="alpha diversity") + 
  stat_boxplot(geom ='errorbar', width=0.2) 

p_richness_RS + stat_compare_means(comparisons = my_comparisons_RS, label = "p.signif" ) + theme_bw() + 
  theme(axis.text.x = element_text(size = 8)) + ggtitle("Richness") + 
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), panel.border = element_blank(),
        panel.background = element_blank(), plot.title = element_text(hjust = 0.5))

p_shannon.e_RS <-  ggplot(alpha_RS_shannon.e , aes(x=Responder, y=value, fill= Responder)) + geom_boxplot(color="black") + 
  labs(x="Responder", y="alpha diversity") + 
  stat_boxplot(geom ='errorbar', width = 0.2)  
p_shannon.e_RS + stat_compare_means(comparisons = my_comparisons_RS, label ="p.signif") + theme_bw() + ggtitle("Shannon Effective") + 
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), panel.border = element_blank(),
        panel.background = element_blank(), plot.title = element_text(hjust = 0.5)) + 
  theme(axis.text.x = element_text(size = 8))

######################################################## GENOTYPE 5WK ##############################################################################


compare_means(value ~ Genotype, data = alpha_gen__5_richness, method = "kruskal.test")
compare_means(value ~ Genotype, data = alpha_gen__5_shannon.e, method = "kruskal.test")

my_comparisons_g5 <- list( c("cre", "wt"))


p_richness_g5 <- ggplot(alpha_gen__5_richness , aes(x=Genotype, y=value, fill= Genotype )) + geom_boxplot(color="black") +  
  labs(x="Genotype", y="alpha diversity") + 
  stat_boxplot(geom ='errorbar', width=0.2) 

p_richness_g5 + stat_compare_means(comparisons = my_comparisons_g5, label = "p.signif" ) + theme_bw() + 
  theme(axis.text.x = element_text(size = 8)) + ggtitle("Richness") + 
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), panel.border = element_blank(),
        panel.background = element_blank(), plot.title = element_text(hjust = 0.5))

p_shannon.e_g5 <-  ggplot(alpha_gen__5_shannon.e , aes(x=Genotype, y=value, fill= Genotype)) + geom_boxplot(color="black") + 
  labs(x="Genotype", y="alpha diversity") + 
  stat_boxplot(geom ='errorbar', width = 0.2)  
p_shannon.e_g5 + stat_compare_means(comparisons = my_comparisons_g5, label ="p.signif") + theme_bw() + ggtitle("Shannon Effective") + 
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), panel.border = element_blank(),
        panel.background = element_blank(), plot.title = element_text(hjust = 0.5)) + 
  theme(axis.text.x = element_text(size = 8))
