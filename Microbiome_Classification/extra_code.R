categorical_variables <- otu_table_scaled_labels[col_name]
cols_factor <- as.factor(otu_table_scaled_labels[[col_name]]) 
levels(cols_factor) <- 1:length(levels(cols_factor))
cols <- as.numeric(cols_factor)
cols <- as.factor(cols)
print(cols)
otu_table_scaled_labels[[col_name]] <- cols
mapping <- data.frame(categorical_variables, otu_table_scaled_labels$Phenotype)
colnames(mapping) <- c("Original", "Encoded")

# function to de-encode categorical variables 
deencoder <- function(df,mapping,col_1, col_2){
  if(missing(col_2)) {
    df[,col_1] <- mapping$Original[match(df[,col_1], mapping$Encoded)]
    return(df)
  } else {
    df[,col_1] <- mapping$Original[match(df[,col_1], mapping$Encoded)]
    df[,col_2] <- mapping$Original[match(df[,col_2], mapping$Encoded)]
    return(df)
  }
}

training = otu_scaled_labels[trainIndex,]
test = otu_scaled_labels[-trainIndex,]

# set X and y - could also be function with above
X_train <- training[,1:(ncol(training)-1)] 
y_train <- training[ , ncol(training)]
X_test <- test[,1:(ncol(test)-1)] 
# for comparing model predictions against actual categories 
actual <- as.factor(test[ , ncol(test)])

mapping_test <- subset(mapping, rownames(mapping) %in% samples)