##################################################################################################################################################
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

# Please select cross-validation method: 
# 0 = k-fold Cross-validation (default) -
# 1 = Repeated k-fold Cross Validation -
# 2 = Leave-one-out cross validation -

cv <- 0       #<--- CHANGE ACCORDINGLY !!!   

# Please give the column where the categorical variable is found 

col_name <- "Phenotype"        #<--- CHANGE ACCORDINGLY !!!   


#############################                           NO CHANGES REQUIRED BELOW THIS LINE                         #############################        

##################################################################################################################################################
###############################################                   Main Script                      ###############################################  
##################################################################################################################################################

###############################################       Load all required libraries       ##########################################################

# Check if required packages are already installed, and install if missing
packages <-c("caret", "ROCR", "dplyr", "xgboost", "pROC") 

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
##################################################################################################################################################
# read data 

otu <- read.table(input_otu_table, sep="\t", header=T, row.names=1, stringsAsFactors=TRUE, comment.char="", check.names=FALSE)
mapping <- read.table(mapping_file, sep="\t", header=T, row.names=1, stringsAsFactors=TRUE, comment.char="", check.names=FALSE)


# scale pre-preprocessed training data and merge phenotype column from metadata

otu_table_scaled <- scale(otu, center = TRUE, scale = TRUE)

otu_table_scaled_labels <- data.frame(t(otu_table_scaled))  
otu_table_scaled_labels[col_name] <- mapping[rownames(otu_table_scaled_labels), col_name]

# convert category to continous variable 
cols <- as.factor(otu_table_scaled_labels[[col_name]]) # make phenotype dynamic, also needs to be some way of informing user which is which
levels(cols) <- 1:length(levels(cols))
cols <- as.numeric(cols)
cols <- as.factor(cols)
print(cols)
otu_table_scaled_labels[[col_name]] <- cols
#otu_table_scaled_labels$Phenotype <- as.integer(as.factor(otu_table_scaled_labels$Phenotype))


# set random seed to 42 
set.seed(42)


#split into training and test 
trainIndex = createDataPartition(otu_table_scaled_labels[[col_name]], 
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

mtryStart <- floor(sqrt(ncol(X_train)))
mtryexpand <- seq(from = mtryStart-10, to= mtryStart+10, by=2) 
tunegrid <- expand.grid(.mtry=mtryexpand)
RF_cv <- train(X_train, y_train, method="rf", ntree=500 , 
                     tuneGrid=tunegrid, trControl=fit_ctrl)
importance <- varImp(RF_cv)
predictions <- predict(RF_cv, newdata = X_test)
prob <- predict(RF_cv, newdata= X_test, type="prob")
samples <- row.names(X_test)
pred_df <- data.frame(samples, actual, predictions) 
print(pred_df)
pred_df$Correct <- pred_df$actual == pred_df$predictions
result <- confusionMatrix(predictions, actual)
metrics <- data.frame(cbind(t(result$positive),t(result$byClass), t(result$overall)))
importance <- importance$importance
otu_names = cbind(OTU=row.names(importance), importance) # clean this code up and make sure it works in all cases (drop index also)
importance_sorted <- otu_names[order(-otu_names$Overall), , drop=FALSE]
pdf("feature_importance_top10.pdf", width = 10, height = 5)
importance_plot <- barplot(importance_sorted[1:10, "Overall"], names.arg=importance_sorted[1:10, "OTU"], 
                           ylab="Variable Importance", las=2, ylim=c(0,100), col = "darkblue",
                           main="Feature Importance (Top 10)")
importance_plot
dev.off()
true_classes <- as.data.frame(test[ , ncol(test)])
model_predictions <- as.data.frame(pred_df$predictions)
write.table(importance, file="importance.tab", sep="\t")
write.table(pred_df, file = "random_forest_predictions.tab", sep="\t", row.names = FALSE) # output to folders
write.table(result$table, file = "confusion_matrix.tab", sep="\t", row.names = FALSE)
write.table(metrics, file="metrics.tab", sep="\t", row.names = FALSE)
if (cv == 0 | 1) {
         #rf_pred = prediction(as.numeric(model_predictions, true_classes))
          #rf.perf = performance(rf.pred,"tpr","fpr")
         # plot(rf.perf,main="ROC Curve for Random Forest",col=2,lwd=2)
         # abline(a=0,b=1,lwd=2,lty=2,col="gray")  
} else if (cv == 2) {
      #  print("Sorry LOOCV support is not yet available")
}









