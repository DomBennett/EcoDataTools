# No copyright
# Author: Dominic John Bennett
# Contains a series of functions for data manipulation of ecological data in R
# (So far only one function)


DegToDec <- function(data) {
  # Takes 'coordinate strings' in degrees (D M (S)) and converts to decimals.
  # Expects coordinate string such: ##'##'##'
  #   - where: # denotes any number of any length, ' denotes any non-numeric(s).
  # Function will stop if more than 3 numbers in coordinate string
  #
  # Args:
  #   data: a matrix or dataframe of coordinates (Latitude and Longitude).
  #
  # Returns:
  #   A dataframe of converted coordinates.
  
  # input checking
  if(!is.data.frame(data) & !is.matrix(data)) {
    stop("Input must be a matrix or dataframe, such:\n
         Latitude     Longitude\n
         51°28\'46.92\"  0°0\'29.47\"")
  }
  if(ncol(data)!= 2) {
    stop("Input must have two columns, such:\n
         Latitude     Longitude\n
         51°28\'46.92\"  0°0\'29.47\"")
  }
  two.numbers <- apply(data, 1:2, function(x) grepl("([^0-9][0-9]|[0-9][^0-9])", x))
  if(!all(two.numbers)) {
    print("Invalid data entries:")
    print(data[!two.numbers])
    stop("Input must be coordinate character strings with at least 2 numbers, e.g.:\n
         Latitude     Longitude\n
         51°28\'46.92\"  0°0\'29.47\"")
  }
  
  # vector inputs
  x <- data[ ,1] #latitude is x
  y <- data[ ,2] #longitude is y
  
  # find start and end of numbers in coordinate strings
  index.x1 <- gregexpr("[0-9\\-\\.]+", x)
  index.x2 <- gregexpr("([0-9][^0-9]|[0-9]$)", x)
  index.y1 <- gregexpr("[0-9\\-\\.]+", y)
  index.y2 <- gregexpr("([0-9][^0-9]|[0-9]$)", y)
  
  # vector outputs
  decimal.x <- rep(NA, nrow(data))
  decimal.y <- rep(NA, nrow(data))
  
  # loop through each row of coordinate string:
  #   - take the first number in string, add the second number/60
  #     and then the third/3600
  #    - stop if more than 3 numbers found.
  error.message1 <- "\nMore than 3 numbers in:\n data["
  error.message2 <- "\nTry checking for non-numeric characters in numbers e.g.:\n
  Here, there is an extra space between the 9 and 2 -- 51°28\'46.9 2\""
  for(i in 1:nrow(data)) {
    decimal.x[i] <- as.numeric(substr(x[i], index.x1[[i]][1], index.x2[[i]][1]))
    for(j in 2:length(index.x1[[i]])) {
      if(j > 3){stop(paste0(error.message1, 1, ",", i,"] -- ",x[i],error.message2))}
      division <- ifelse(j == 2, 60, 3600)
      next.x <- as.numeric(substr(x[i], index.x1[[i]][j], index.x2[[i]][j]))
      decimal.x[i] <- decimal.x[i] + (next.x/division)
    }
    decimal.y[i] <- as.numeric(substr(y[i], index.y1[[i]][1], index.y2[[i]][1]))
    for(j in 2:length(index.y1[[i]])) {
      if(j > 3){stop(paste0(error.message1, 2, ",", i,"] -- ",y[i],error.message2))}
      division <- ifelse(j == 2, 60, 3600)
      next.y <- as.numeric(substr(y[i], index.y1[[i]][j], index.y2[[i]][j]))
      decimal.y[i] <- decimal.y[i] + (next.y/division)
    }
  }
  
  # control for S and W
  decimal.x <- decimal.x * ifelse(grepl("(S|s)", x), -1, 1)
  decimal.y <- decimal.y * ifelse(grepl("(W|w)", y), -1, 1)
  
  # generate output
  decimals <- data.frame(Latitude = decimal.x, Longitude = decimal.y)
  return(decimals)
}