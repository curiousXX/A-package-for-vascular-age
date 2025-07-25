#' Title Data Processing for the model construction
#' Description This function processes the input data frame by checking the data types of its features and converting them to numeric if they meet certain conditions. It also counts the number of missing values in each row and column, and returns a list containing continuous features, discrete features, and NA counts.
#'
#' @param Dat A dataframe containing the variable to be processed.
#' @param pd A positive integer that is used to determine the minimum number of unique values required for a variable to be considered continuous. The default value is 5.
#' @param na_threshold A numeric value between 0 and 1 that represents the threshold for the proportion of missing values allowed in a variable. The default value is 0.5.
#'
#' @returns A dataframe with the same structure as the input data frame, but with the data types of the specified features converted to numeric.
#' @export
#' @import lubridate
#' @importFrom stats na.omit prcomp quantile sd varimax
#' @examples
#'
#' data(test_data)
#'
#' # Process the test_data with default parameters
#' Dat <- data_process(test_data)
#'
data_process <- function(Dat, pd = 5, na_threshold = 0.5) {
  # Check if the loaded data is a data frame
  if (!inherits(Dat, "data.frame")) {
    message("Input validation: Dat is not a data frame, attempting to convert to a data frame")
    original_type <- class(Dat)[1]
    # Check if the input is a matrix, table, or data.table
    allowed_types <- c("matrix", "table", "data.table")
    if (!original_type %in% allowed_types) {
      stop(paste("Unsupported input type:", original_type))
    }
  }  else {
    message("Input validation: Dat is a data frame")
  }

  # Attempt to convert the input to a data frame and check for errors
  # Use tryCatch to handle potential errors during conversion
  tryCatch({
    Dat <- as.data.frame(Dat)
  }, error = function(e) {
    # Enhanced error handling
    error_message <- paste("Data type conversion failed:",
                           "Original type:", original_type,
                           "Error:", e$message)
    stop(error_message)
  })
  # Check if the loaded data has at least one column
  if (ncol(Dat) == 0) {
    stop("Dat must have at least one column.")
  }
  # Check if the loaded data has at least one row
  if (nrow(Dat) == 0) {
    stop("Dat must have at least one row.")
  }
  # Enhanced colnames check: Check for duplicate column names
  if (any(duplicated(colnames(Dat)))) {
    dup_cols <- colnames(Dat)[duplicated(colnames(Dat))]
    warning(paste("Duplicate column names detected:", paste(dup_cols, collapse = ", ")))
  }
  # Check if pd is a positive integer
  if (!is.numeric(pd) || pd <= 0 || pd != round(pd)) {
    stop("pd must be a positive integer.")
  }
  # Check if na_threshold is between 0 and 1
  if (!is.numeric(na_threshold) || na_threshold < 0 || na_threshold > 1) {
    stop("na_threshold must be between 0 and 1.")
  }
  # The detect_and_convert_date() function is used to identify and convert date columns in the dataset
  # It checks each column in the dataset to see if it can be converted to a date format
  # If a column can be converted, it is converted to the Date type and stored in the original dataset
  detect_and_convert_date <- function(Dat) {
    # Define a function to clean strings by replacing certain characters and trimming whitespace
    clean_strings <- function(x) {
      x <- gsub("[/\\]", "-", x)  # Uniform separators
      x <- gsub("^\\s+|\\s+$", "", x)  # Remove the first and last spaces
      x
    }
    date_cols <- c()
    for (col in colnames(Dat)) {
      # The clean_strings() is used to preprocess the strings before attempting to convert them to dates
      cleaned_col <- clean_strings(Dat[[col]])
      # 扩展日期格式支持
      orders <- c("ymd", "mdy", "dmy", "ydm", "myd", "dym",
                  "ymd HMS", "mdy HMS", "dmy HMS")
      temp <- suppressWarnings(lubridate::parse_date_time(
        cleaned_col,
        orders = orders,
        exact = FALSE
      ))
      # Add the secondary verification of date validity
      if (all(!is.na(temp)) &&
          length(unique(temp)) > 1 &&  # Exclude all the same dates
          !all(lubridate::year(temp) == 1900)) {  # Exclude invalid default dates
        Dat[[col]] <- temp
        date_cols <- c(date_cols, col)
      }
    }
    if (length(date_cols) > 0) {
      message(paste0("Date columns detected and converted: ", paste(date_cols, collapse = ", ")))
    } else {
      message("No date columns detected.")
    }
    return(Dat)
  }
  # Call the detect_and_convert_date() function to convert date columns in the dataset
  Dat <- detect_and_convert_date(Dat)
  # The C.check() function is made to filter dataset by checking several specific conditions related to completeness and variability
  C.check <- function(x, pd, na_threshold) {
    # Check if x is a numeric vector or can be coerced to numeric
    x_numeric <- tryCatch({
      as.numeric(na.omit(as.character(x)))
    }, warning = function(w) {
      return(NULL)
    }, error = function(e) {
      return(NULL)
    })
    # If x_numeric is NULL, return 0
    if (is.null(x_numeric)) return(0)
    # Check if the data is constant, i.e., all values are the same. .Machine$double.eps^0.5 is a small number that represents the smallest difference between two distinct numbers in R.
    # If the standard deviation is less than this value, we consider it constant.
    is_constant <- sd(x_numeric, na.rm = TRUE) < .Machine$double.eps^0.5
    # Check if the data is outlier, i.e., if any value is more than 5 standard deviations away from the mean.
    outlier_test <- any(abs(scale(x_numeric)) > 5, na.rm = TRUE)
    # Three specific conditions have to be TRUE to return 1, if one of them is FALSE then return 0
    # 1. Number of NAs in x is less than na_threshold of its total length
    # 2. Number of unique values is more than pd
    # 3. The data is not constant
    if (length(unique(x_numeric)) > pd &&
        !is_constant &&
        sum(is.na(x_numeric)) < (length(x_numeric) * na_threshold)) {
      return(1)
    }
    0
  }
  # This code snippet checks whether each column in the data frame `Dat` meets certain conditions
  # defined by the function `C.check`. The `apply` function is used to apply `C.check` column-wise.
  # Parameters passed to `C.check` include:
  # - `x`: A column from the data frame `Dat`.
  # - `pd`: A parameter used within the `C.check` function.
  # - `na_threshold`: A threshold value for handling missing data.
  # The result of the `apply` function is converted into a data frame `is_contin`.
  # The resulting data frame has a single column named 'result'.
  is_contin <- as.data.frame(apply(Dat, 2, function(x) C.check(x, pd, na_threshold)))
  colnames(is_contin) <- 'result'
  # The fsite() function is used to find the indices of specific values in the column names of a matrix or data frame.
  fsite <- function(v, m) {
    # Extracts the column names of m and assigns them to the 'site.name' variable
    site.name = colnames(m)
    # Make an empty vector named site to store the indices of the matched columns
    site = c()
    # Make a loop with number of iteration equal to length of v
    # If there's an index in site.name match with one element of v, assign it to x
    # Then x is appended and accumulated in the site vector
    for(i in 1:length(v)) {
      x = which(site.name == v[i])
      site = c(site, x)
    }
    return(site)
  }
  # Use the fsite() function to find columns identified as continuous variables
  # Continuous variables are defined for each TRUE (1) value in column "is_contin"
  site <- fsite(row.names(is_contin)[is_contin$result == 1], Dat)
  # If there are continuous variables, they will be put in vector contin_fea
  ifelse(length(site) > 0,
         # If there are continuous variables, assign their column names to contin_fea
         contin_fea <- colnames(Dat)[site],
         # If there are no continuous variables, assign an empty string to contin_fea
         contin_fea <- '')
  # Use the fsite() function to extract columns identified as non-continuous variables
  # Non-continuous variables are defined for each FALSE (0) value in column "is_contin"
  site <- fsite(row.names(is_contin)[is_contin$result == 0], Dat)
  # Extract columns identified as non-continuous variables and put them into a new data frame 'discrete_Dat'
  # If there are non-continuous variables, they will be put in a new data frame 'discrete_Dat'
  if (length(site) > 0) {
    # Extract columns identified as non-continuous variables and put them into a new data frame 'discrete_Dat'
    discrete_Dat <- Dat[, site, drop = F]
  }
  # If there are no non-continuous variables, assign blank spaces to both variables 'discrete_fea_two' and 'discrete_fea_more'
  if (length(site) == 0) {
    discrete_fea_two <- discrete_fea_more <- ''
    # Make a list consisting of variables 'contin_fea', 'discrete_fea_two', and 'discrete_fea_more' to show the Continuous Feature and Discrete Feature
    result <- list(contin_fea,
                   discrete_fea_two,
                   discrete_fea_more)
    # Name each element in the list
    names(result) <- c('Contin_feature',
                       'Discrete_feature(level<=2)',
                       'Discrete_feature(level>2)')
    # Return the result list
    return(result)
  }
  # The str2() function is mainly used to separate discrete and multi-class variables.
  str2 <- function(x){
    # Convert the data frame into character, then remove the NA values
    x <- na.omit(as.character(x))
    # Convert all elements in x into factor, so we can treat unique values as levels
    # Get all the unique levels of the factor, then count the number of unique levels
    # Check if the length (amount) of unique levels is two or fewer
    # If TRUE, return 1. If FALSE, return 0.
    if(length(levels(factor(x))) <= 2)
    {
      return(1)
    }
    return(0)
  }
  # Apply function str2() to the columns of data frame discrete_Dat (parameter 1 means row, 2 means column)
  # Then convert the output into a data frame and assign it to the is_two variable
  is_two <- as.data.frame(apply(discrete_Dat, 2, function(x) str2(x)))
  # Rename the is_two column to "result"
  colnames(is_two) <- 'result'
  # Use the fsite() function to find columns identified as discrete variables with 2 or fewer levels
  site <- fsite(row.names(is_two)[is_two$result == 1], discrete_Dat)
  # If there are discrete variables with 2 or fewer levels, they will be put in vector discrete_fea_two
  ifelse(length(site) > 0,
         # If there are discrete variables with 2 or fewer levels, assign their column names to discrete_fea_two
         discrete_fea_two <- colnames(discrete_Dat)[site],
         # If there are no discrete variables with 2 or fewer levels, assign an empty string to discrete_fea_two
         discrete_fea_two <- '')
  # Use the fsite() function to find columns identified as discrete variables with more than 2 levels
  site <- fsite(row.names(is_two)[is_two$result == 0], discrete_Dat)
  # If there are discrete variables with more than 2 levels, they will be put in vector discrete_fea_more
  ifelse(length(site) > 0,
         # If there are discrete variables with more than 2 levels, assign their column names to discrete_fea_more
         discrete_fea_more <- colnames(discrete_Dat)[site],
         # If there are no discrete variables with more than 2 levels, assign an empty string to discrete_fea_more
         discrete_fea_more <- '')
  # Put all the variables in the list and name them
  # Return the list of continuous, discrete (level <= 2), and discrete (level > 2) features
  result <- list(contin_fea,
                 discrete_fea_two,
                 discrete_fea_more)
  names(result) <- c('Contin_feature',
                     'Discrete_feature(level<=2)',
                     'Discrete_feature(level>2)')
  # The na_count() function is used to count the number of NA (missing) values for each row and column
  # Define a function named na_count that takes a dataset 'Dat' and a boolean 'result' as arguments
  na_count <- function(Dat,result = TRUE){
    # Calculate and print the overall proportion of NA values in the entire dataset
    message(paste0('NA ratio : ',sum(is.na(Dat))/(ncol(Dat)*nrow(Dat))))
    # -----Count NA Values for Each Row-----
    # Use function apply() to count the number of NA using function sum(is.na(x))
    # margin=1 indicates that the function should be applied to rows
    df_row <- apply(Dat,1,function(x){sum(is.na(x))})
    # Convert the obtained data into data frame
    df_row <- as.data.frame(df_row)
    # Name the column as 'count', then convert it into numeric format
    colnames(df_row) <- 'count'
    df_row$count <- as.numeric(as.character(df_row$count))
    # Calculate the proportion of missing values in each row to three decimal places
    df_row$perc <- round((df_row$count)/ncol(Dat),3)
    # Sort the data frame in descending order based on the number of NA values
    df_row <- df_row[order(df_row$count,decreasing = TRUE),]
    # -----Count NA Values for Each Column-----
    # Use function apply() to count the number of NA using function sum(is.na(x))
    # margin=2, then function count NA values for each column
    df_col <- apply(Dat,2,function(x){sum(is.na(x))})
    # Convert the obtained data into data frame
    df_col <- as.data.frame(df_col)
    # Name the column as 'count', then convert it into numeric format
    colnames(df_col) <- 'count'
    df_col$count <- as.numeric(as.character(df_col$count))
    # Calculate the proportion of missing values in each column to three decimal places
    df_col$perc <- round((df_col$count)/nrow(Dat),3)
    # Sort the data frame in descending order based on the number of NA values
    df_col <- df_col[order(df_col$count,decreasing = TRUE),]
    # Make a list concsist of both NA values in each rows and columns data frames, then assinged into 'result_output' variable
    result_output <- list(df_row,df_col)
    # Give each of the list a name, then return the 'result_output' variable as function output
    names(result_output) <- c('row','col')
    return(result_output)
  }
  # Call the na_count() function to count the number of NA values in the dataset 'Dat'
  # Count NA using the na_count() function, and stored in 'na_result' variable
  na_result <- na_count(Dat)
  # Merge both result from data types determination and NA amount counting in 'result' list
  result <- list(contin_fea,
                 discrete_fea_two,
                 discrete_fea_more,
                 na_result)
  # Name each element in the list
  names(result) <- c('Contin_feature',
                     'Discrete_feature(level<=2)',
                     'Discrete_feature(level>2)',
                     'NA_count')
  # The bat_as.change() function is used to convert the data types of the specified features in the dataset
  bat_as.change <- function(name.v,Dat){
    # Put input dataset into temporary variable
    temp_df <- Dat
    # Find features in input variable 'temp_df' that match the data from input 'name.v' variable
    # The result is stored in 'site.v' variable using fsite() function
    site.v <- fsite(name.v,
                    temp_df)
    # Make a loop with number of iteration equal to length of 'site.v'
    # If there's an index in site.name match with one element of 'site.v', convert it into numeric value
    # Then the converted value is appended back to variable 'temp_df'
    for(i in 1:length(site.v)){
      temp_df[,site.v[i]] <- as.numeric(as.character(temp_df[,site.v[i]]))
    }
    return(temp_df)
  }
  # Extract the continuous variables from 'fea_result' variable
  v.contin <- result$Contin_feature
  new_df <- bat_as.change(v.contin,Dat)
  # logic transformation for continuous variables
  log_transform <- function(new_Dat, v.contin) {
    transform_log <- function(x) {
      # Try to apply log transformation and handle errors
      tryCatch({
        if (all(x > 0, na.rm = TRUE)) {
          log(x)
        } else {
          # Adjust the values to be positive before applying log transformation
          adjusted <- x - min(x, na.rm = TRUE) + 1e-6
          log(adjusted)
        }
      }, error = function(e) {
        warning(paste("Log transformation failed:", e$message))
        x
      })
    }
    for (col in v.contin) {
      message(paste("Applying log transformation to column:", col))
      new_col_name <- paste0("log_", col)
      new_Dat[[new_col_name]] <- transform_log(new_Dat[[col]])
      # Remove the new column if all values are NA
      if (all(is.na(new_Dat[[new_col_name]]))) {
        new_Dat[[new_col_name]] <- NULL
        warning(paste("All values are NA in column:", new_col_name, "Column removed."))
      }
      # Check for infinite values and replace them with NA
      if (any(is.infinite(new_Dat[[new_col_name]]))) {
        new_Dat[[new_col_name]][is.infinite(new_Dat[[new_col_name]])] <- NA
        warning(paste("Infinite values detected in column:", new_col_name, "Replaced with NA."))
      }
    }
    return(new_Dat)
  }
  # Apply log transformation to continuous variables
  if (length(v.contin) > 0) {
    new_df <- log_transform(new_df, v.contin)
  } else {
    message("No continuous variables to apply log transformation.")
  }
  # Return the final result list containing continuous features, discrete features, and NA counts
  return(new_df)
}
