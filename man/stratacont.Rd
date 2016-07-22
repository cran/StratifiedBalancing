\name{stratacont}
\alias{stratacont}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Stratacont()
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This is the primary stratification function for continuous outcome variables.It locates naturally occurring strata in the data, weights them and returns an improved regression coefficient for the treatment variable which is adjusted for both confounding and interaction terms. Returns a matrix containig naturally occurring strata.
}
\usage{
stratacont(Treatment,Outcome,Matrix)
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
This is the primary stratification function for continuous outcome variables.It locates naturally occurring strata in the data, weights them and returns an improved regression coefficient for the treatment variable which is adjusted for both confounding and interaction terms.
}

\examples{
## We will first begin by simulating data in 5 covariates and a continuous outcome
## with significant interaction terms and correlations amongst covariates (to simulate an
## experiment with a strongly non-linear underlying model).
## First, we will create a matrix with the input variables. The inout variables will all be
## categorical variables.
m=matrix(nrow=5000,ncol=6)
for ( i in 1:ncol(m)){
m[,i]=rbinom(5000,1,0.5)
}


## Next, we will simulate the output variable and include interaction terms
for(i in 1:nrow(m)){
a=(2*m[i,5] + 0.5*m[i,1] + 4*m[i,2] + 2.3*m[i,3] + 5*m[i,4] +
2.3*m[i,3]*m[i,2] +3.5*m[i,1]*m[i,2] + 2.1*m[i,1]*m[i,3] +
5*m[i,1]*m[i,2]*m[i,3] + 6*m[i,1]*m[i,4] +3*m[i,2]*m[i,4] +
2*m[i,3]*m[i,4] + 3.4*m[i,1]*m[i,2]*m[i,3]*m[i,4] +
5*m[i,1]*m[i,2]*m[i,4] + 4*m[i,2]*m[i,3]*m[i,4])
m[i,6]=rnorm(1,a,1)
}
## We are interested in determining the coefficient of covariate 5 which is 2.
## Tmost straightforward
## way of doing this is to use simple linear regression as follows
m=as.data.frame(m)
k=lm(m[,6]~.,data=m[,(1:5)])
## The value of the coefficient of variable 5 found by the regression can be retrieved using
k$coeff[6]

## We can now use the stratacont() function to find a more accurate estimation of the coefficient
g=stratacont(5,6,m)
## Note that as the model includes more covariates, the accuracy of the stratification
## techniques is far superior.

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
