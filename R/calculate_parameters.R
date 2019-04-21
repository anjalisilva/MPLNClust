calculate_parameters<-function(g,y){
  d<-ncol(y)
  mu_para<-d*g
  sigma_para<-(d*((d+1)/2))*g
  pi_para<-g-1 # because if you have g-1 parameters, you can do 1-these to get the last one
  paratotal<-mu_para+sigma_para+pi_para # total parameters are
  return(paratotal)
}

