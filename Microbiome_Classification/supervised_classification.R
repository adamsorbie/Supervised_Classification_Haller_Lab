##############################################################################################################################################
# This script performs supervised classification of a normalised OTU table. 
# Version: 0.9 
# 2018-02-05  Author: Adam Sorbie 
#
#  Please set the directory of the script as the working folder (e.g D:/studyname/NGS-Data/Rhea/beta-diversity/)
#' Note: the path is denoted by forward slash "/"
setwd("C:/Users/PhD/Supervised_Classification/Microbiome_Classification")  #<--- CHANGE ACCORDINGLY !!!

# Enter name of OTU table file: 
input_otu_table <- "test_otu.tab"        #<--- CHANGE ACCORDINGLY !!!

# Enter name of mapping file: 
mapping_file <- "test_map.tab"
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
packages <-c("caret", "randomForest", "ROCR") 

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

otu <- read.table(input_otu_table, sep="\t", header=T, row.names=1, stringsAsFactors=FALSE, comment.char="", check.names=FALSE)
mapping <- read.table(mapping_file, sep="\t", header=T, row.names=1, stringsAsFactors=FALSE, comment.char="", check.names=FALSE)

# scale pre-preprocessed training data (normalised relative abundance filtered 1% abundance in at least one sample) and merge phenotype column from metadata

otu_table_scaled <- scale(otu, center = TRUE, scale = TRUE)

otu_table_scaled_labels <- data.frame(t(otu_table_scaled))  
otu_table_scaled_labels[col_name] <- mapping[rownames(otu_table_scaled_labels), col_name] 

# set random seed to 42 
set.seed(42)

# set X and y 
X <- otu_table_scaled_labels[,1:(ncol(otu_table_scaled_labels)-1)] 
y <- otu_table_scaled_labels[ , ncol(otu_table_scaled_labels)]


if (model == 0) {
  RF_phenotype_classify <- randomForest( x=X , y=y , ntree=500, mtry = c(1:13), importance=TRUE, proximities=TRUE  )
} else if (model == 1) {
  print("Unfinished")
} else if (model == 2) {
  print("Unfinished")
} else if (model == 3) {
  print("Unfinished")
} else { 
  print("fail")}
  




