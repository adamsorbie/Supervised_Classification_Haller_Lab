##################################################################################################################################################
# This script performs supervised classification of a normalised OTU table. 
# Version: 0.6.0
# 2020-01-06  Author: Adam Sorbie 
#
#  Please set the directory of the script as the working folder (e.g C:/studyname/NGS-Data/Microbiome_classification)
#' Note: the path is denoted by forward slash "/"
#' Alternatively if you open this script from the working directory comment out the code below with a preceding hashtag (#)
setwd("C:/Users/PhD/Supervised_Classification/Microbiome_Classification")  #<--- CHANGE ACCORDINGLY !!!

# Enter name of OTU table file: 
input_otu_table <- "merged_otu.tab"        #<--- CHANGE ACCORDINGLY !!!

# Enter name of mapping file: 
mapping_file <- "merged_map-madeupmulticlass.tab"         #<--- CHANGE ACCORDINGLY !!!

# Please select cross-validation method: 
# 0 = k-fold Cross-validation (default) -
# 1 = Repeated k-fold Cross Validation -
# 2 = Leave-one-out cross validation -

cv <- 0    #<--- CHANGE ACCORDINGLY !!!   

# Please give the column where the categorical variable is found 

col_name <- "Phenotype"        #<--- CHANGE ACCORDINGLY !!!   


#############################                           NO CHANGES REQUIRED BELOW THIS LINE                         #############################        

##################################################################################################################################################
###############################################                   Main Script                      ###############################################  
##################################################################################################################################################

###############################################       Load all required libraries       ##########################################################

# could update this to pacman or at least add if else for bioconductor 
# Check if required packages are already installed, and install if missing
packages <-c("caret", "dplyr", "xgboost", "pROC", "ggplot2", "compositions") 

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
#################################################################################################################################################
# functions - places where you do things twice could possibly be re-written with applies

# might not need this
# defaults <- list(model = "Random Forest", 
#                  train_size = 0.7,
#                  )

#' to-do
#' 1. comment functions properly
#' 2. check all function outputs
#' 3. error handling
#' 4. plotting function 
#' 5. include abundance in plotting function and return plot with subplot
#' 6. write output function 
#' 7. return any other useful vars
#' 8. add other model(s)
#' 9. ROC curve plotting function 
#' 10. LOOCV support (needs for loop)
#' 11. ML reading make sure everything is ok 
#' 12. full function documentation 
#' 13. source functions
#' 14. write as package with accompanying script


preprocess <- function(otu, scale = "clr", mapping, class_col) {
  # scale otu using centred log ratio transform 
  scaled_otu <- clr(otu)
  # add class column (y) to otu table
  otu_scaled_labels <- data.frame(t(scaled_otu))
  otu_scaled_labels["Class"] <- mapping[rownames(otu_scaled_labels), class_col]
  
  return(otu_scaled_labels)
}

check_dim <- function(df) {
  # if sample size is left than 50 print warning
  if (dim(otu)[2] < 50) {
    print("Warning: small sample size detected, classifications and feature importances may be less accurate/useful")
  }
}

get_cv <- function(cv) {
  # select cross-validation method 
  if (cv == 0) {
    fit_ctrl <- trainControl(method = "cv", number = 10)  
  } else if (cv == 1) {
    fit_ctrl <- trainControl(method = "repeatedcv", number = 10, repeats = 3)
  } else if (cv == 2) {
    fit_ctrl <- trainControl(method = "LOOCV")
  } else {
    # if no cv method selected return error
    print("Error please enter a valid selection:
        0 = k-fold cross validation
        1 = repeated k-fold cross validation
        2 = leave-one-out cross validation")
  }
  return(fit_ctrl)
}

train_test_split <- function(otu_scaled_labels, class_col, partition) {
  # generate indices for splitting data into test and training
  trainIndex <- createDataPartition(otu_scaled_labels[[class_col]], p=partition,
                                    list = F, times = 1)
  # 
  training <- otu_scaled_labels[trainIndex, ]
  ncol_training <- ncol(training)
  test <- otu_scaled_labels[-trainIndex,]
  ncol_test <- ncol(test)
  
  actual <- as.factor(test[ , ncol_test])
  
  # generate X and y - probably a much cleaner way to write this 
  X_train <- training[,1:(ncol_training-1)]
  y_train <- training[ , ncol_training]
  X_test <- test[, 1:(ncol_test -1)]
  
  # return list
  return_list <- list("X_train" = X_train, "y_train" = y_train, "X_test" = X_test, "Actual" = actual, 
                      "training_full" = training, "test_full" = test)
  return(return_list)
}

rf_tune <- function(X, mtry_step) {
  
  mtryStart <- floor(sqrt(ncol(X)))
  
  mtryexpand <- seq(from = mtryStart-mtry_step, to= mtryStart +mtry_step, by = 2)
  tunegrid <- expand.grid(.mtry=mtryexpand)
  return(tunegrid)
}

# xgboost_tune <- function(X, )

build_model <- function(X, y, method, ...){
  model <- train(X, y, method = method, trControl = fit_ctrl, ...)
  return(model)
}



get_importance_df <- function(model, topx) {
  imp <- varImp(model)
  imp <- imp$importance 
  
  if (missing(topx)) {
    imp_sort <- imp[order(imp$Overall, decreasing = T), , drop=F]
  }
  else{
    imp_sort <- imp[order(imp$Overall, decreasing = T), , drop=F] %>% 
      top_n(topx)
  }
  
  return(imp_sort)
}

pred_stats <- function(names, actual, model_predictions) {
  df <- data.frame(names, actual, model_predictions)
  df$Correct <- df$actual == df$model_predictions
  return(df)
}

return_model_output <- function(model, X_test, output_prob = NULL, actual, ...) {
  importance <- get_importance_df(model, ...)
  predictions <- predict(model, newdata=X_test)
  if (output_prob == TRUE){
    prob <- predict(model, newdata=X_test, type="prob")
  }
  confusion_matrix <- confusionMatrix(predictions, actual)
  pred_statistics <- pred_stats(row.names(X_test), actual = actual, model_predictions = predictions)
  
  return_list <- list("Variable_Importance" = importance,
                      "Model_predictions" = predictions,
                      "Confusion_matrix" = confusion_matrix, 
                      "Prediction_accuracy" = pred_statistics)
  if (exists("prob")) {
    return_list <- c(return_list, list("Probabilities" = prob))
  }
  return(return_list)
}

# 
# if (cv == 0 | 1) {
#   if (length(unique(categorical_variables[[col_name]])) == 2) {
#     roc_rf <- roc(actual, prob$`1`)
#     pdf("roc_curve.pdf")
#     roc_plot <- plot(roc_rf, col = "blue")
#     roc_plot 
#     dev.off() 
#   }   else {
#     roc_rf <- multiclass.roc(actual, prob$`1`)
#     pdf("roc_curve.pdf")
#     roc_plot <- plot(roc_rf$rocs[[3]], col = "blue")
#     roc_plot
#     dev.off() 
#   } 
# } else if (cv == 2) {
#   print("Sorry LOOCV support is not yet available") 
# }


##################################################################################################################################################
# read data 

otu <- read.table(input_otu_table, sep="\t", header=T, row.names=1, stringsAsFactors=TRUE, 
                  comment.char="", check.names=FALSE)
mapping <- read.table(mapping_file, sep="\t", header=T, row.names=1, stringsAsFactors=TRUE, 
                      comment.char="", check.names=FALSE)

# check dimensions of otu table and print warning message if sample size small 
check_dim(otu)
stopifnot(exists("otu"), exists("mapping"))
# set random seed to 42 
set.seed(42)

# scale pre-preprocessed training data and merge phenotype column from metadata - should maybe think about using different method of scaling clr?
otu_scaled_labels <- preprocess(otu, scale = "clr", mapping, col_name)

#split into training and test 
train_test <- train_test_split(otu_scaled_labels, class_col = "Class", partition = 0.7)

# cross validation method  - re-do as function input cv and return fit_ctrl
fit_ctrl <- get_cv(cv)

# 
tunegrid <- rf_tune(train_test$X_train, mtry_step = 10)
RF_cv <- train(train_test$X_train, train_test$y_train, trControl=fit_ctrl, method="rf", ntree=500 , 
                     tuneGrid=tunegrid)

# return model importance and predictions 
model_out <- return_model_output(RF_cv, X_test = train_test$X_test, output_prob = F, actual = train_test$Actual)



metrics <- data.frame(cbind(t(result$positive),t(result$byClass), t(result$overall)))

pdf("feature_importance_top10.pdf", width = 10, height = 5)

# write ggplot plotting function to plot var importance with abundance subplot 

# ggplot(importance_10, aes(x= "OTU", y="Overall") ) + geom_bar(stat="identity", fill = "indianred") + 
       # coord_flip() + ggtitle("Feature Importance (Top 10")
importance_plot <- barplot(importance_sorted[1:10, "Overall"], names.arg=importance_sorted[1:10, "OTU"], 
                           ylab="Variable Importance", las=2, ylim=c(0,100), col = "darkblue",
                           main="Feature Importance (Top 10)")
importance_plot
dev.off()
true_classes <- as.data.frame(test[ , ncol(test)])
model_predictions <- as.data.frame(pred_df$predictions)

# writing files could be wrapped into function 
write.table(importance, file="importance.tab", sep="\t")
write.table(pred_df, file = "random_forest_predictions.tab", sep="\t", row.names = FALSE) # output to folders
write.table(result$table, file = "confusion_matrix.tab", sep="\t", row.names = FALSE)
# write.table(metrics, file="metrics.tab", sep="\t", row.names = FALSE)

# this can also be function 











