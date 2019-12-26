# Parameter calculation
calcParameters <- function(numberG, dimensionality) {
  
  muPara <- dimensionality * numberG
  
  sigmaPara <- (dimensionality * ((dimensionality + 1) / 2)) * numberG
  # because if you have numberG-1 parameters, 
  # you can do 1-these to get the last one
  
  piPara <- numberG - 1 
  
  # total parameters
  paraTotal <- muPara + sigmaPara + piPara
  
  return(paraTotal)
  # Developed by Anjali Silva
}
