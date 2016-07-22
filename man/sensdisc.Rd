\name{sensdisc}
\alias{sensdisc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
sensdisc()

}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function performs sensitivity analysis on discrete outcome variables
by dropping varaibles one at a time, in-order to determine which variables can be removed without altering the results of stratified covariate balancing.
}
\usage{
sensdisc(Treatment,Outcome,Matrix)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Treatment}{
%%     ~~Describe \code{x} here~~
Column number of variable to be used as treatment.
}
\item{Outcome}{
%%     ~~Describe \code{x} here~~
Column number of variable to be used as outcome.
}
\item{Matrix}{
%%     ~~Describe \code{x} here~~
Name of matrix or data.frame where data is stored.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
This function performs sensitivity analysis by dropping varaibles one at a time, in-order to determine which variables can be removed without altering the results of stratified covariate balancing.
}
\examples{
## In this example we will generate a matrix with a large number of
## covariates and a small number of observations. No model will be
## built into the data, our goal here is to demonstrate how sensitivty
## analysis would be performed.

## Firstly a matrix with 10 columns and 1000 observations will be created
m=matrix(nrow=1000,ncol=10)
for(i in 1:ncol(m)){
m[,i]=rbinom(1000,1,0.5)
}

## Next we will run function sensdisc() on the data.
g=sensdisc(9,10,m)

}

\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
