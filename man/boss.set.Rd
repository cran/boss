\name{boss.set}
\Rdversion{1.1}
\alias{boss.set}
%\alias{genCor}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Create an object to be passed to boss.fit
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Takes as arguments everything in a GWAS model that does not explicitely involve genotype and pre-computes as much as possible
}
\usage{
boss.set(formula, E.name = NULL, family = gaussian(), id = NULL, 
	corstr = "independence", type = "glm", method = "chol", data, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
	A formula, of the kind used by the 'type' specified, that does not include genotype.
}
  \item{E.name}{
  For interaction studies, the name (as a character) of the variable with which to perform the interaction
}
  \item{family}{
  A family object, only gaussian() and binomial() supported
}
  \item{id}{
  A vector of unique id's for subjects in the sample, to be matched with id's associated with genotype.  For GEE models, also used to generate clusters.
}
  \item{corstr}{
  For GEE models only, one of "independence", "ar1" and "exchangeable", specifying the working correlation structure
}
  \item{type}{
  One of "glm", "gee" or "lmm"
}
  \item{method}{
   One of "chol" or "swap", default is "chol".  "swap" cannot be used with robust standard errors or interaction studies, but is substantially faster.
}
  \item{data}{
  An optional data frame in which the variables in `formula' are located.
}
  \item{\dots}{
  Arguments to be passed to the `glm()', 'lmer()' or 'geese()' used internally in boss.set. The argument `sattdf=TRUE' can be used for
}
}

\value{
An object to be used in `boss.set()'
}
\references{
Voorman, A.; Rice, K. and Lumley, T. \emph{Fast Computation for Genome Wide Association Studies using Boosted One-Step Statistics}. Bioinformatics 2012; doi: 10.1093/bioinformatics/bts291 

Satterwaite F. \emph{An Approximate Distribution of Estimates of Variance Components}. Biometrics Bulletin 1946; 2(6) 110-114.

Pan W, Wall M. \emph{Small-sample adjustments in using the sandwich variance estimator in generalized estimating equations}. Statist Med 2002; 21:1429-1441.
}
\author{
Arend Voorman, Colleen Sitlani
}
\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{boss.fit}}
\code{\link{boss.set}}
\code{\link{glm}}
\code{\link{geese}}
\code{\link{lmer}}
}
\examples{
#generate data:
n <- 500
nsnps <- 1000
snpMat <- replicate(nsnps,rbinom(n,2,.01))
colnames(snpMat) <- 1e6+(1:nsnps)
rownames(snpMat) <- 1:n

#generate 3 observations per person with random intercept and covariates:
id <- as.numeric(gl(n,3))
y <- rnorm(3*n) + rnorm(n,sd=.3)[id]
data <- data.frame(y=y, x1=rnorm(3*n,mean=y), x2 =rnorm(3*n), 
	site=sample(6,n,replace=TRUE)[id], id=id )

# set up a GEE with independence working correlation structure, using "chol" method:
init1 <- boss.set(y~x1+x2+factor(site), id=data$id, type = "gee", 
	method = "chol", corstr= "independence", data=data)

# set up a linear mixed model (lmm) with random intercept using "swap" method:
init2 <- boss.set(y~x1+x2+factor(site) + (1|id), type = "lmm", method = "swap", data=data)

#set up an interaction with logistic regression, which must use the "chol" method:
y <- rbinom(n,1,0.5) 
data.bin <- data.frame(y=y, x1=rnorm(n,mean=y), x2 =rnorm(n), 
	e = rnorm(n), site=sample(6,n,replace=TRUE), id=1:n)

init3 <- boss.set(y~x1+x2+factor(site)+e, E.name = "e", id=id, 
	family=binomial(), type="glm", data = data.bin)

##fit the model (at one locus):
g <- snpMat[,1]

###GEE
#one step approximation, using satterwaite degrees of freedom.
re <- boss.fit(g[id], init1, robust=TRUE, sattdf=TRUE)
re

#p-value with df approximation:
2*pt(abs(re$beta.main/sqrt(re$v.main)),df=re$df, lower=FALSE)
#p-value with normal approximation:
2*pnorm(abs(re$beta.main/sqrt(re$v.main)), lower=FALSE)


#full iteration
boss.fit(g[id], init1, robust=TRUE, thresh=1)


##LMM
#one step approximation:
boss.fit(g[id], init2)

#full iteration (uses lmer, since "swap" method does not permit further iterations currently)
mod <- lmer(update(init2$formula,.~.+g[id]),data=data)
list(beta.main = getME(mod,"beta")[9], v.main=vcov(mod)[9,9],
     chi2 = getME(mod,"beta")[9]^2/vcov(mod)[9,9])

##Logistic Reg.
boss.fit(g, init3)

#full iteration
boss.fit(g, init3, thresh = 1)

#for swap method, matrix arguments can be used to simultaneously fit many models:
results <- boss.fit(t(snpMat[,id]),init2)


}
%\keyword{ ~htest }
%\keyword{ ~regression }
