data{
  int<lower=1> d; // Dimension of theta
  int<lower=0> N; //Sample size
  int y[N,d]; //Array of Y
  vector[d] mu;
  matrix[d,d] Sigma;
  vector[d] normfactors;  
}
parameters{
  matrix[N,d] theta;
}
model{ //Change values for the priors as appropriate
  for (n in 1:N){
    theta[n,]~multi_normal(mu,Sigma);
  }
  for (n in 1:N){
    for (k in 1:d){
      real z;
      z=exp(normfactors[k]+theta[n,k]);
      y[n,k]~poisson(z);
    }
  }
}


