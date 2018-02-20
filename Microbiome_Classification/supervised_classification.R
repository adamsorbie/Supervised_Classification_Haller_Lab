##############################################################################################################################################
# This script performs supervised classification of a normalised OTU table. 
# Version: 0.9 
# 2018-02-05  Author: Adam Sorbie 
#
#  Please set the directory of the script as the working folder (e.g C:/studyname/NGS-Data/Microbiome_classification)
#' Note: the path is denoted by forward slash "/"
setwd("C:/Users/PhD/Supervised_Classification/Microbiome_Classification")  #<--- CHANGE ACCORDINGLY !!!

# Enter name of OTU table file: 
input_otu_table <- "merged_otu.tab"        #<--- CHANGE ACCORDINGLY !!!

# Enter name of mapping file: 
mapping_file <- "merged_map.tab"         #<--- CHANGE ACCORDINGLY !!!
# Please select model.
# 0 = Random-Forest Model (default) - 
# 1 = Support Vector Machine - 
# 2 = eXtreme Gradient Boosting -   
 
model <- 0       #<--- CHANGE ACCORDINGLY !!!

# Please select cross-validation method: 
# 0 = k-fold Cross-validation (default) -
# 1 = Repeated k-fold Cross Validation -
# 2 = Leave-one-out cross validation -

cv <- 1       #<--- CHANGE ACCORDINGLY !!!   

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

# scale pre-preprocessed training data and merge phenotype column from metadata

otu_table_scaled <- scale(otu, center = TRUE, scale = TRUE)

otu_table_scaled_labels <- data.frame(t(otu_table_scaled))  
otu_table_scaled_labels[col_name] <- mapping[rownames(otu_table_scaled_labels), col_name] 

# set random seed to 42 
set.seed(42)


#split into training and test 
trainIndex = createDataPartition(otu_table_scaled_labels$Phenotype, 
                                 p=0.7, list=FALSE,times=1)
training = otu_table_scaled_labels[trainIndex,]
test = otu_table_scaled_labels[-trainIndex,]

# set X and y 
X_train <- training[,1:(ncol(training)-1)] 
y_train <- training[ , ncol(training)]
X_test <- test[,1:(ncol(test)-1)] 
# for comparing model predictions against actual categories 
actual <- as.character(test[ , ncol(test)])

# cross validation method  
if (cv == 0) {
   fit_ctrl <- trainControl(method = "cv", number = 10)  
} else if (cv == 1) {
    fit_ctrl <- trainControl(method = "repeatedcv", number = 10, repeats = 3)
} else if (cv == 2) {
    fit_ctrl <- trainControl(method = "LOOCV")
} else {
  print("Error please enter a valid selection:
        0 = k-fold cross validation
        1 = repeated k-fold cross validation
        2 = leave-one-out cross validation")
}


if (model == 0) {
    mtryStart <- floor(sqrt(ncol(X_train)))
    mtryexpand <- seq(from = mtryStart-10, to= mtryStart+10, by=2) 
    tunegrid <- expand.grid(.mtry=mtryexpand)
    RF_cv <- train(X_train, y_train, method="rf", ntree=501 , 
                   tuneGrid=tunegrid, trControl=fit_ctrl)
    importance <- varImp(RF_cv)
    predictions <- predict(RF_cv, newdata = X_test)
    samples <- row.names(X_test)
    pred_df <- data.frame(samples, actual, predictions) 
    print(pred_df)
    pred_df$Correct <- pred_df$actual == pred_df$predictions
    result <- confusionMatrix(predictions, actual)
    metrics <- data.frame(cbind(t(result$positive),t(result$byClass), t(result$overall)))
    importance <- importance$importance
    otu_names = cbind(OTU=row.names(importance), importance)
    importance_sorted <- arrange(importance, otu_names  , desc(Overall)  )

    
    importance_plot
    write.table(importance, file="importance.tab", sep="\t")
    write.table(pred_df, file = "random_forest_predictions.tab", sep="\t", row.names = FALSE) # output to folders
    write.table(result$table, file = "confusion_matrix.tab", sep="\t", row.names = FALSE)
    write.table(metrics, file="metrics.tab", sep="\t", row.names = FALSE)
} else if (model == 1) {
  svm_cv <- train(X_train, y_train, method="svmLinear", trControl= fit_ctrl ) #tunegrid 
  predictions <- predict(svm_cv, newdata = X_test)
  samples <- row.names(X_test)
  pred_df <- data.frame(samples, actual, predictions) 
  print(pred_df)
  pred_df$Correct <- pred_df$actual == pred_df$predictions
  result <- confusionMatrix(predictions, actual)
  print(svm_cv)
} else if (model == 2) {
  print("Sorry, this model is not yet available, please choose another")
} else { 
  print("Error, please enter a valid selection: 
        0 = Random Forest
        1 = SVM
        2 = eXtreme gradient boosting")}
  




