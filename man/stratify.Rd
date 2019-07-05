\name{stratify}
\alias{stratify}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Function Stratify
%%  ~~function to do ... ~~
}
\description{
This is the main function which perfroms Stratified Covariate Balancing. It also enables the user to opt for stratifying only the parents in the Markov Blanket of the Treatment variable along with discretization and use synthetic cases to match the data.
}
\usage{
stratify(Treatment, Outcome, Matrix, Discretize, Synthetic, Ordered , Markov)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Treatment}{
The column number of the Treatment variable
}
  \item{Outcome}{
The column number of the outcome variable
}
  \item{Matrix}{
The name of the data.frame or matrix where the data is stored
}
  \item{Discretize}{
A TRUE/FALSE parameter indicating whether covariates should be discretized or not
}
  \item{Synthetic}{
A TRUE/FALSE parameter indicating whether synthetic matching should be used or not
}
  \item{Ordered}{
A TRUE/FALSE parameter indicating whether the variables provided are ordered or not
}
  \item{Markov}{
A TRUE/FALSE parameter indicating whether Markov Blanket of the treatment should be used for rstratification or not
}
}
\examples{

#create snthetic data
m=matrix(nrow=100,ncol=5,data=0)
for(i in 1:ncol(m)){
m[,i] = rbinom(100,1,0.5)
}
colnames(m)=c("Var1" , "Var2" , "Var3" , "var4" , "Var5")

## Now use Stratify
g=stratify(4,5,m)
}





