##############################################################################################################################################
# This script performs supervised classification of a normalised OTU table. 
# Version: 0.9 
# 2018-02-05  Author: Adam Sorbie 
#
#  Please set the directory of the script as the working folder (e.g D:/studyname/NGS-Data/Rhea/beta-diversity/)
#' Note: the path is denoted by forward slash "/"
setwd("C:/Users/PhD/Supervised_Classification/Microbiome_Classification")  #<--- CHANGE ACCORDINGLY !!!

# Enter name of OTU table file: 
input_otu_table <- "merged_otu.tab"        #<--- CHANGE ACCORDINGLY !!!

# Enter name of mapping file: 
mapping_file <- "merged_map.tab"
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
packages <-c("caret", "randomForest", "ROCR", "plyr", "rfUtilities") 

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

otu <- read.table(input_otu_table, sep="\t", header=T, row.names=1, stringsAsFactors=TRUE, comment.char="", check.names=FALSE)
mapping <- read.table(mapping_file, sep="\t", header=T, row.names=1, stringsAsFactors=TRUE, comment.char="", check.names=FALSE)

# scale pre-preprocessed training data (normalised relative abundance filtered 1% abundance in at least one sample) and merge phenotype column from metadata

otu_table_scaled <- scale(otu, center = TRUE, scale = TRUE)

otu_table_scaled_labels <- data.frame(t(otu_table_scaled))  
otu_table_scaled_labels[col_name] <- mapping[rownames(otu_table_scaled_labels), col_name] 

# set random seed to 42 
set.seed(42)


#split into training and test 
trainIndex = createDataPartition(otu_table_scaled_labels$Phenotype, 
                                 p=0.7, list=FALSE,times=1)
train = otu_table_scaled_labels[trainIndex,]
test = otu_table_scaled_labels[-trainIndex,]

# set X and y 
X_train <- train[,1:(ncol(train)-1)] 
y_train <- train[ , ncol(train)]
X_test <- test[,1:(ncol(test)-1)] 
# for comparing predictions against actual categories 
actual <- test[ , ncol(test)]
# cross validation method  
if (cv == 0) {
    train_ctrl <- trainControl(method = "repeatedcv", number = 10, repeats = 3)  
} else if (cv == 1) {
    train_ctrl <- trainControl(method = "repeatedcv", number = 10, repeats = 3)
} else if (cv == 2) {
    train_ctrl <- trainControl(method = "repeatedcv", number = 10, repeats = 3)
} else {
  
}


if (model == 0) {
    mtryStart <- floor(sqrt(ncol(X))) 
    mtry_test <- tuneRF(X_train, y_train, mtryStart, ntreeTry=500, 
                        stepFactor=2, improve=0.05, trace=TRUE, plot=TRUE)
    mtry <- as.data.frame(mtry_test)
    mtry_sorted <- mtry[order(mtry$OOBError),]
    mtry_use <- mtry_sorted$mtry[1]
    RF_classify <- randomForest(X_train, y_train, ntree=500, 
                                mtry = mtry_use, importance=TRUE, proximities=TRUE  )
    print(RF_classify)
    predictions <- predict(RF_classify, newdata = X_test)
    print(predictions)
} else if (model == 1) {
  print("Sorry, this model is not yet available, please choose another")
} else if (model == 2) {
  print("Sorry, this model is not yet available, please choose another")
} else if (model == 3) {
  print("Sorry, this model is not yet available, please choose another")
} else { 
  print("Error, please enter a valid choice: 
        0 = Random Forest
        1 = SVM
        2 = eXtreme gradient boosting
        3 = Neural network (Multilayer Perceptron")}
  




