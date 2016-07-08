\name{stratadisc}
\alias{stratadisc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Stratadisc()

}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This is the primary stratification function for discrete outcome variables.It locates naturally occurring strata in the data, weights them and returns a common odds ratio for the treatment variable which is adjusted for both confounding and interaction
terms.Returns a matrix containing naturally occurring strata.
}
\usage{
stratadisc(Treatment,Outcome,Matrix)
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
This is the primary stratification function for discrete outcome variables.It locates naturally occurring strata in the data, weights them and returns a common odds ratio for the treatment variable which is adjusted for both confounding and interaction terms.
}

\examples{
## We will first begin by simulating data in 11 covariates and a discrete outcome, with
## significant interaction terms and correlations amongst covariates (to simulate a
## non-randomized experiment with a strongly non-linear underlying model).
## First, we will create a matrix with the input variables
m=matrix(nrow=1000,ncol=12)
for ( i in 1:ncol(m)){
m[,i]=rbinom(1000,1,0.5)
}
## Next, we will create correlations amongst covariates 2,4 and 11
for( i in 1:nrow(m)){
if(m[i,11] == 1)m[i,2]=rbinom(1,1,0.8)
if(m[i,11] == 0)m[i,2]=rbinom(1,1,0.2)
if(m[i,11] == 1)m[i,4]=rbinom(1,1,0.8)
if(m[i,11] == 0)m[i,4]=rbinom(1,1,0.2)
}
## Next, we will simulate the output variable and include interaction terms
for(i in 1:nrow(m)){
a=exp(2*m[i,11] + 0.5*m[i,1] - 4*m[i,2] + 2.3*m[i,3] +
0.8*m[i,4] -0.7*m[i,5] - 4*m[i,6] + 3.6*m[i,7] +
1.2*m[i,8] - 11*m[i,9] - 2.1*m[i,10] + 2.3*m[i,3]*m[i,4]
-3.5*m[i,5]*m[i,6]*m[i,7] + 8*m[i,1]*m[i,2]*m[i,9] +
2.1*m[i,2]*m[i,6]*m[i,8]) / (1 + exp(2*m[i,11] + 0.5*m[i,1]
- 4*m[i,2] + 2.3*m[i,3] + 0.8*m[i,4] -0.7*m[i,5] - 4*m[i,6]
+ 3.6*m[i,7] + 1.2*m[i,8] - 11*m[i,9] - 2.1*m[i,10] +
2.3*m[i,3]*m[i,4] -3.5*m[i,5]*m[i,6]*m[i,7] +
8*m[i,1]*m[i,2]*m[i,9] + 2.1*m[i,2]*m[i,6]*m[i,8]))
m[i,12]=rbinom(1,1,a)
}
## We are interested in determining the coefficient of covariate 11 which is 2.
##The most straightforward way of doing this
##is to use logistic regression as follows
m=as.data.frame(m)
k=glm(m[,12]~.,data=m[,(1:11)],family=binomial)
## The value of the coefficient of variable 11 found by the
##logistic regression can be retrieved using
k$coeff[12]

## We can now use the stratadisc() function to find a
##more accurate estimation of the coefficient
g=stratadisc(11,12,m)
## We need to take the log() of the first
##number returned "Odds Ratio of Impact Of Treatment On Outcome"
## and the estimated value of the parameter should be more accurate

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
