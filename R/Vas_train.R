#' Title Model construction based on KDM after standardization and PCA
#'
#' @param Dat A dataframe containing the variable to be processed.
#' @param v.name The biomarkers used in the model you want to construct.
#' @param gender_col The colname of gender in your data frame. The default value is "gender".
#' @param male A character representing male in the gender column. The default value is "male".
#' @param female A character representing female in the gender column. The default value is "female"
#' @param split A numeric value indicating the split method for standardization. split can be 1, 2, or 3. If split is 1, all data will be standardized and the model of all data will constructed. If split is 2, the data will be standardized and model will be constructed after splitting by gender. If split is 3, you can obtain all the models. The default value is 2.
#' @param age_col The name of age column in the data frame. The default value is "age".
#'
#' @returns A list containing the following elements: the standardization parameters of the biomarkers used in the train model you want to construct, which will be used in the test data. Moreover, the model using the KDM method will be constructed.
#' @export
#' @importFrom stats na.omit prcomp quantile sd varimax
#' @import bioage
#'
#' @examples
#'
#' data(test_data)
#'
#' # You can start by organizing the data with data_process function.
#' # Calculate the missing data rate and identify continuous variables for log transformation
#' Dat <- data_process(test_data)
#'
#' # Example usage of Vas_calcu function
#' # v.name is the colname of biomarkers used in the model you want to construct.
#' # If the column specified by v.name contains a missing value, the function will stop.
#' try_result <-Vas_train(Dat, v.name = c("PWV", "ABI","MAP"), gender_col = "gender", male = "male",female = "female", split = 2, age_col = "age")
#'
#' # Extract the standardization parameters of the biomarkers used in the train model you want to construct
#' try_result$std_params
#'
#' # Extract the KDM model constructed using the data you input
#' # The parameter of this model are used to calculate the biological age of the input data
#' train_model <- try_result$train_model
#'
#' # Extract the male vascular age
#' vascular_age <- try_result$age_pred$male$data$bioage
#'
#' # Extract the female vascular age
#' vascular_age <- try_result$age_pred$female$data$bioage
#'
#' # Extract the difference between vascular age and actual age
#' vascular_age_advanced <- try_result$age_pred$male$data$baaccel_diff
#'
#'
Vas_train <- function(Dat, v.name = c("MAP", "PWV", "ABI"), gender_col = "gender", male = "male", female = "female", split = 2, age_col = "age") {
  # Check if the loaded data is a data frame
  if (!inherits(Dat, "data.frame")) {
    message("Input validation: Dat is not a data frame, attempting to convert to a data frame")
    original_type <- class(Dat)[1]
    tryCatch({
      Dat <- as.data.frame(Dat)
      message(paste("Successfully converted", original_type, "type to data frame"))
    }, error = function(e) {
      stop("Data type conversion failed:", e$message)
    })
  } else {
    message("Input validation: Dat is already a data frame")
  }
  # Check if the loaded data has at least 10 rows
  # If the number of rows is less than 10, print a warning message
  if (nrow(Dat) < 10) {
    warning("Small dataset warning: Sample size is less than 10, results may be unstable")
  } else {
    message("Input validation: Sample size is sufficient")
  }
  # Check the number of the specified biomarkers
  # If the number of biomarkers is less than 2, stop the function and print a message
  # If the number of biomarkers is greater than or equal to 2, print a message indicating that the number of biomarkers is valid
  if (length(v.name) < 2) {
    stop("The number of biomarkers must be at least 2.")
  } else {
    message("The number of biomarkers is valid.")
  }
  # Check for outliers in the specified biomarkers
  for (var in v.name) {
    quartiles <- quantile(Dat[[var]], probs = c(0.25, 0.75), na.rm = TRUE)
    IQR <- IQR(Dat[[var]], na.rm = TRUE)
    lower_bound <- quartiles[1] - 1.5 * IQR
    upper_bound <- quartiles[2] + 1.5 * IQR
    outliers <- Dat[[var]][Dat[[var]] < lower_bound | Dat[[var]] > upper_bound]
    if (length(outliers) > 0) {
      message(paste("Outlier check: The following outliers were found in", var, ":", paste(outliers, collapse = ", ")))
    }
  }
  # Check if the loaded data has the required columns
  required_columns <- c(v.name, gender_col, age_col)
  # Check if the required_columns are present in the data frame
  missing_columns <- setdiff(required_columns, colnames(Dat))
  # If any required_columns are missing, stop the function and print the missing columns
  # If all required_columns are present, print a message indicating that all required_columns are present
  if (length(missing_columns) > 0) {
    stop(paste("The following biomarkers are missing from the data frame:", paste(missing_columns, collapse = ", ")))
  } else {
    message("All required_columns are present in the data frame.")
  }
  # Validate that all biomarker columns are numeric
  # If any biomarker columns are not numeric, stop the function and print a message
  # If all biomarker columns are numeric, print a message indicating that the biomarker columns are valid
  non_numeric_biomarkers <- v.name[!sapply(Dat[, v.name, drop = FALSE], is.numeric)]
  if (length(non_numeric_biomarkers) > 0) {
    stop("The following biomarkers are not numeric:", paste(non_numeric_biomarkers, collapse = ", "))
  }
  # Check if the specified biomarkers contain missing values
  # If any missing values are found, stop the function and print a message
  # If no missing values are found, print a message indicating that there are no missing values
  na_report <- sapply(Dat[, v.name], function(x) {
    na_count <- sum(is.na(x))
    if(na_count > 0) paste(na_count, "missing values") else "no"
  })
  if (any(grepl("missing values", na_report))) {
    message("Missing value check: The following biomarkers contain missing values:")
    print(data.frame(Variable = names(na_report), Status = na_report))
    stop("Please check the missing values in the specified biomarkers.")
  } else {
    message("Missing value check: No missing values in the specified biomarkers.")
  }
  # Enhanced gender classification validation
  if (!is.factor(Dat[[gender_col]])) {
    message("Gender column is not factor type, attempting to convert...")
    tryCatch({
      Dat[[gender_col]] <- as.factor(Dat[[gender_col]])
    }, error = function(e) {
      stop("Failed to convert gender column to factor: ", e$message)
    })
  }
  gender_levels <- unique(Dat[[gender_col]])
  if (!all(c(male, female) %in% gender_levels)) {
    message("The detected gender categories in the data are -> ", paste(gender_levels, collapse = ", "))
  }
  if (length(gender_levels) > 2) {
    message("Detected more than 2 gender categories: ", paste(gender_levels, collapse = ", "))
    gender_counts <- table(Dat[[gender_col]])
    print(gender_counts)
  }
  # Check if the age column is numeric and within a reasonable range (0-120)
  if (!is.numeric(Dat[[age_col]])) {
    message("Age column is not numeric, attempting to convert...")
    tryCatch({
      Dat[[age_col]] <- as.numeric(Dat[[age_col]])
    }, error = function(e) {
      stop("Failed to convert age column to numeric: ", e$message)
    })
  }
  age_values <- Dat[[age_col]]
  if (any(age_values < 0 | age_values > 120)) {
    invalid_ages <- age_values[age_values < 0 | age_values > 120]
    message("The number of abnormal age values: ", length(invalid_ages))
    message("The abnormal age values are: ", paste(invalid_ages, collapse = ", "))
    warning("Unusual age values detected (<0 or >120). Please check data quality.")
  }
  # Check if split is numeric 1， 2, or 3
  if (!is.numeric(split)) {
    stop("split must be numeric")
  }
  if (split %% 1 != 0) {
    split <- as.integer(split)
    message("split is not an integer, converting to integer")
  }
  if (!split %in% 1:3) {
    stop("split must be 1, 2, or 3")
  }
  # The split parameter is automatically adjusted based on the number of gender factors
  if (length(gender_levels) == 1){
    split <- 2
    message("Gender has only one level; Only one level result is output")
  }
  # The auto_std() function is used to standardize the data
  auto_std <- function(Dat, v.name, gender_col, male, female, split){
    # Initialize an empty list to store standardization tables for each gender
    tab_std_list <- list()
    standardized_data <- list()
    # Calculate the standardization parameters (mean and standard deviation) for each variable in v.name
    # The calc_std_params() function is used to calculate the mean and standard deviation of the specified variables
    # The function takes a data frame and a vector of variable names as input and returns a data frame with the mean and standard deviation for each variable
    calc_std_params <- function(data_subset, v.name) {
      params <- data.frame(
        term = v.name,
        u = sapply(data_subset[, v.name], mean, na.rm = TRUE),
        sd = sapply(data_subset[, v.name], sd, na.rm = TRUE)
      )
      return(params)
    }
    # The apply_standardization() function is used to standardize the data using the mean and standard deviation
    # The function takes a data frame and a data frame of standardization parameters as input and returns the standardized data frame
    # The function loops through each variable name in the standardization parameters and applies the standardization formula
    # If the standard deviation is not equal to 0, standardize the data by subtracting the mean and dividing by the standard deviation
    # If the standard deviation is equal to 0, standardize the data by subtracting the mean only
    # The function returns the standardized data frame
    apply_standardization <- function(data_subset, params) {
      for (i in seq_along(params$term)) {
        term <- params$term[i]
        u <- params$u[i]
        sd_val <- params$sd[i]
        if (sd_val != 0) {
          data_subset[, term] <- (data_subset[, term] - u) / sd_val
        } else {
          data_subset[, term] <- data_subset[, term] - u
        }
      }
      return(data_subset)
    }
    # Check if the split is 1 or 3, all data is used for standardization
    if (split %in% c(1, 3)) {
      # Standardize the data using all data
      Dat_all <- Dat
      # Using the calc_std_params() function to calculate the mean and standard deviation of the specified variables
      tab_std_all <- calc_std_params(Dat, v.name)
      # Using the apply_standardization() function to standardize the data using the mean and standard deviation
      Dat_all <- apply_standardization(Dat_all, tab_std_all)
      # Store the standardization parameters and standardized data in the list
      tab_std_list$all <- tab_std_all
      standardized_data$all <- Dat_all
    }
    # Check if the split is 2 or 3, feamle and male data are used for standardization
    if (split %in% c(2, 3)){
      # Check if there is male data in the gender column
      if (any(Dat[, gender_col] == male)) {
        # Standardize the male data
        Dat_male <- Dat[Dat[[gender_col]] == male, ]
        # Using the calc_std_params() function to calculate the mean and standard deviation of the specified variables
        tab_std_male <- calc_std_params(Dat_male, v.name)
        # Using the apply_standardization() function to standardize the data using the mean and standard deviation
        Dat_male <- apply_standardization(Dat_male, tab_std_male)
        # Store the standardization parameters and standardized data in the list
        tab_std_list$male <- tab_std_male
      } else {
        Dat_male <- NULL
        tab_std_list$male <- NULL
      }
      # Check if there is female data in the gender column
      if (any(Dat[, gender_col] == female)) {
        Dat_female <- Dat[Dat[,gender_col] == female, ]
        # Using the calc_std_params() function to calculate the mean and standard deviation of the specified variables
        # Using the apply_standardization() function to standardize the data using the mean and standard deviation
        tab_std_female <- calc_std_params(Dat_female, v.name)
        Dat_female <- apply_standardization(Dat_female, tab_std_female)
        # Store the standardization parameters and standardized data in the list
        tab_std_list$female <- tab_std_female
      } else {
        Dat_female <- NULL
        tab_std_list$female <- NULL
      }
      # Rbind the data
      if (!is.null(Dat_male) && !is.null(Dat_female)) {
        standardized_data$gender <- rbind(Dat_male, Dat_female)
      } else if (!is.null(Dat_male)) {
        standardized_data$gender <- Dat_male
      } else if (!is.null(Dat_female)) {
        standardized_data$gender <- Dat_female
      }
    }
    # Return the standardization tables and the processed data frame
    return(list(params = tab_std_list, Data = standardized_data))
  }
  # The 'run_pca()' function is used to perform PCA analysis on the data
  run_pca <- function(data, gender_col, age_col, v.name_pc) {
    # Check if the data frame is empty, and if so, print a message and return NULL
    if (is.null(data) || nrow(data) == 0) {
      message("The data frame is empty. PCA cannot be performed.")
      return(NULL)
    }
    # The specification of the PCA analysis is as follows:
    # center=F -> the data is not centered, centered by subtracting with the mean
    # scale=F -> the data is not scaled (usually, scaled by dividing with the standard deviation)
    # The PCA analysis is performed on the specified variables in v.name
    pca <- prcomp(data[, v.name], center = FALSE, scale = FALSE)
    # Calculate the standard deviation of the principal components and divide it by the sum of the standard deviations. This will give the proportion of variance explained by each principal component.
    # The eigenvector is the rotation matrix of the principal components
    pca$sdev/sum(pca$sdev)
    eigenvector <- pca$rotation
    # Calculate the loadings of the principal components
    # The loadings are the eigenvectors multiplied by the standard deviations and stored in the 'rawLoadings' variable
    rawLoadings <- pca$rotation[,1:ncomp] %*% diag(pca$sdev, ncomp, ncomp)
    # The original loadings are rotated using the varimax() function to obtain simpler and more interpretable factor loadings
    # The varimax() function returns the rotated loading matrix, which can be extracted through $loadings.
    rotatedLoadings <- varimax(rawLoadings)$loadings
    # Convert the principal components of the training set to a data frame
    dat_pca <- as.data.frame(pca$x)
    # Merge the 'Age' column from the 'data_std_train' data frame to the 'dat_pca' data frame
    dat_pca <- cbind(dat_pca, data[, c(gender_col, age_col)])
    # Construct the bioage model using the principal components
    # The 'bioage::kdm_calc()' function is used to calculate the bioage model using the principal components
    # The function takes the data frame, age variable, and biomarker names as input and returns the bioage model
    # The function returns the bioage model
    model <- bioage::kdm_calc(dat_pca, agevar = age_col, biomarkers = v.name_pc)
    age_pred <- bioage::kdm_calc(dat_pca, agevar = age_col, biomarkers = v.name_pc,fit = model$fit)
    # Return the PCA results and the bioage model
    return(list(pca = pca, model = model, eigenvector = eigenvector,dat_pca = dat_pca, age_pred = age_pred))
  }
  # Call the auto_std() function to standardize the data
  list_std <- auto_std(Dat, v.name, gender_col, male, female, split)
  # df_std is the standardization parameters
  df_std <- list_std[[1]]
  # The standardized train set is stored in the 'data_std_train' variable
  data_std_train <- list_std[[2]]
  # The number of components used in the basic model is equal to the number of basic model indicators
  ncomp <- length(v.name)
  # v.name_pc is the name of the principal components in the model and is stored in the 'v.name_pc' variable
  v.name_pc <- paste0("PC",c(1:length(v.name)))
  # Extract the basic model indicators from the standardized training sets
  train_model <- list()
  eigenvector <- list()
  age_pred <- list()
  eta <- list()
  # Train a combined model using all participants when split is 1 or 3
  if (split %in% c(1, 3)) {
    # Extract the standardized training dataset containing all participants
    t.train <- data_std_train$all

    # Select relevant columns: gender, age, and feature variables for model training
    # These indicators will be used to build the baseline model
    t.train_all <- t.train[,c(gender_col, age_col, v.name)]

    # Run Principal Component Analysis (PCA) on the combined dataset
    train_model$all <- run_pca(t.train_all, gender_col, age_col, v.name_pc)

    # Extract eigenvectors (principal component directions) from the trained model
    eigenvector$all <- train_model$all$eigenvector

    # Extract age prediction results from the trained model
    age_pred$all <- train_model$all$age_pred

    # Calculate eta effect size: (actual age variance - brain age error variance) / actual age variance
    # Higher eta values indicate better model performance
    eta$all <- (age_pred$all$data$agevar - age_pred$all$data$ba.e) / age_pred$all$data$agevar
  }

  # Train separate gender-specific models when split is 2 or 3
  if (split %in% c(2, 3)) {
    # Extract the gender-stratified standardized training dataset
    t.train <- data_std_train$gender

    # Select relevant columns: gender, age, and feature variables
    t.train <- t.train[,c(gender_col, age_col, v.name)]

    # Split the dataset by gender to create male and female subsets
    t.train_male <- t.train[t.train[,gender_col] == male, ]
    t.train_female <- t.train[t.train[,gender_col] == female, ]

    # Train model for male participants if male data is available
    if (nrow(t.train_male) > 0) {
      # Run PCA analysis on male-only dataset
      train_model$male <- run_pca(t.train_male, gender_col, age_col, v.name_pc)

      # Extract eigenvectors from the male-specific model
      eigenvector$male <- train_model$male$eigenvector

      # Extract age prediction results from the male-specific model
      age_pred$male <- train_model$male$age_pred

      # Calculate eta effect size for the male-specific model
      eta$male <- (age_pred$male$data$agevar - age_pred$male$data$ba.e) / age_pred$male$data$agevar
    }

    # Train model for female participants if female data is available
    if (nrow(t.train_female) > 0) {
      # Run PCA analysis on female-only dataset
      train_model$female <- run_pca(t.train_female, gender_col, age_col, v.name_pc)

      # Extract eigenvectors from the female-specific model
      eigenvector$female <- train_model$female$eigenvector

      # Extract age prediction results from the female-specific model
      age_pred$female <- train_model$female$age_pred

      # Calculate eta effect size for the female-specific model
      eta$female <- (age_pred$female$data$agevar - age_pred$female$data$ba.e) / age_pred$female$data$agevar
    }
  }
  # Return a comprehensive list containing all trained models and related data structures
  # train_std: standardized training data used for model development
  # train_model: trained PCA models (all/male/female specific models)
  # data_std_train: standardized training datasets organized by group
  # eigenvector: principal component directions/loadings for each trained model
  # age_pred: age prediction results and performance metrics for each model
  # eta: effect size measurements for evaluating model prediction accuracy
  return(list(train_std = df_std, train_model = train_model, data_std_train = data_std_train, eigenvector = eigenvector, age_pred = age_pred, eta = eta))
}
