##############################################################################################################################################
# This script performs supervised classification of a normalised OTU table. 
# Version: 0.9 
# 2018-02-05  Author: Adam Sorbie 
#
#  Please set the directory of the script as the working folder (e.g D:/studyname/NGS-Data/Rhea/beta-diversity/)
#' Note: the path is denoted by forward slash "/"
setwd("C:/Users/PhD/Supervised_Classification/Microbiome_Classification")  #<--- CHANGE ACCORDINGLY !!!

# Enter name of OTU table file: 
input_otu_table <- "OTUs_w_metadata-Norm.tab"        #<--- CHANGE ACCORDINGLY !!!
  
# Please select model.
# 0 = Random-Forest Model (default) - 
# 1 = Support Vector Machine - 
# 2 = eXtreme Gradient Boosting - 
# 3 = Neural Network (Multi-layer Perceptron) - 
 
model <- 0        #<--- CHANGE ACCORDINGLY !!!

# Please select cross-validation method: 
# 0 = k-fold Cross-validation (default) -
# 1 = Repeated k-fold Cross Validation -
# 2 = Leave-one-out cross validation -

cv <- 0        #<--- CHANGE ACCORDINGLY !!!   

# Please give the column where the categorical variable is found 

col_name <- "Phenotype"        #<--- CHANGE ACCORDINGLY !!!   


  
######                  NO CHANGES REQUIRED BELOW THIS LINE                 ######

##################################################################################
######                             Main Script                              ###### 
##################################################################################

###################       Load all required libraries     ########################

# Check if required packages are already installed, and install if missing
packages <-c("caret", "RandomForest", "ROCR") 

# Function to check whether the package is installed
InsPack <- function(pack)
{
  if ((pack %in% installed.packages()) == FALSE) {
    install.packages(pack,repos = "http://cloud.r-project.org/")
  } 
}

# Applying the installation on the list of packages
lapply(packages, InsPack)

# Make the libraries
lib <- lapply(packages, require, character.only = TRUE)

# Check if it was possible to install all required libraries
flag <- all(as.logical(lib))

# read data 

df <- read.table(input_otu_table)

# scale pre-preprocessed training data (normalised relative abundance filtered 1% abundance in at least one sample) and merge phenotype column from metadata

otu_table_scaled <- scale(df, center = TRUE, scale = TRUE)

otu_table_scaled_labels <- data.frame(t(otu_table_scaled))  
otu_table_scaled_labels$col_name <- metadata[rownames(otu_table_scaled_labels), col_name] 



if (model == 0) {
  
} else if (model == 1) {
  
} else if (model = 2) {
  
} else if (model == 3) {

} else { 
}
  

