#takes degrees Degrees Minutes and Seconds(optional) and converts it to decimals
#expects coordinate string such: ##'##'##', where: # denotes any number of any length, ' denotes any non-numeric(s)
#deals with N, S, E and W

deg2dec <- function(data, warnings = T){
  #unambiguous names
  X <- data[,1] #latitude is X
  Y <- data[,2] #longitude is Y
  #Find numbers in coordinate strings
  indexX1 <- gregexpr("[0-9\\-\\.]+", X)
  indexX2 <- gregexpr("([0-9][^0-9]|[0-9]$)", X)
  indexY1 <- gregexpr("[0-9\\-\\.]+", Y)
  indexY2 <- gregexpr("([0-9][^0-9]|[0-9]$)", Y)
  #empty vectors to collect decimal outpus
  decimalX <- rep(NA, nrow(data))
  decimalY <- rep(NA, nrow(data))
  #loop through each row of coordinate string:
    #take the first number in string, add the second number/60 and then the third/3600
    #Stop if more than 3 numbers found.
  error.message1 <- "\nMore than 3 numbers in:\n data["
  error.message2 <- "\nTry checking for non-numeric characters between coordinates e.g. 3°45'5 4"
  for(i in 1:nrow(data)){
    decimalX[i] <- as.numeric(substr(X[i], indexX1[[i]][1], indexX2[[i]][1]))
    for(j in 2:length(indexX1[[i]])){
      if(j > 3 & warnings){stop(paste0(error.message1, 1, ",", i,"] -- ",X[i],error.message2))}
      division <- ifelse(j == 2, 60, 3600)
      nextX <- as.numeric(substr(X[i], indexX1[[i]][j], indexX2[[i]][j]))
      decimalX[i] <- decimalX[i] + (nextX/division)
    }
    decimalY[i] <- as.numeric(substr(Y[i], indexY1[[i]][1], indexY2[[i]][1]))
    for(j in 2:length(indexY1[[i]])){
      if(j > 3 & warnings){stop(paste0(error.message1, 2, ",", i,"] -- ",Y[i],error.message2))}
      division <- ifelse(j == 2, 60, 3600)
      nextY <- as.numeric(substr(Y[i], indexY1[[i]][j], indexY2[[i]][j]))
      decimalY[i] <- decimalY[i] + (nextY/division)
    }
  }
  #control for S and W
  decimalX <- decimalX * ifelse(grepl("S", X), -1, 1)
  decimalY <- decimalY * ifelse(grepl("W", Y), -1, 1)
  #generate output
  decimals <- data.frame(Latitude = decimalX, Longitude = decimalY)
  return(decimals)
}





##HI DOM
