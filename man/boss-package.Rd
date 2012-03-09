\name{boss-package}
\Rdversion{1.1}
\alias{boss-package}
\alias{boss}
\docType{package}
\title{
~~Fast one-step approximations for GLM, GEE and Mixed Models~~
}
\description{
BOSS uses parameter estimates obtained without genotype to boost standard one-step approximations, and precomputes as much as possible without genotype (using boss.set) to minimize effort required. 
}
\details{
boss.set & pre-computes as much as possible without genotype \cr
boss.fit &  fits a model with genotype, using an object from boss.set \cr
boss.ncdf & is a wrapper for interface with the 'ncdf' package \cr
}
\author{
Arend Voorman

Maintainer: Arend Voorman <voorma@u.washington.ed>
}
\references{
\emph{Fast Computation for Genome Wide Association Studies using Boosted One-Step Statistics}
Arend Voorman; Ken Rice; Thomas Lumley
Bioinformatics 2012; doi: 10.1093/bioinformatics/bts291

\emph{Performance of generalized estimating equations in practical situations}
Lipsitz, S.R. and Fitzmaurice, G.M. and Orav, E.J. and Laird, N.M.
Biometrics, 1994, pp270-278
}


\keyword{ package }
\examples{
n <- 500
nsnps <- 1000
snpMat <- replicate(nsnps,rbinom(n,2,.2))
colnames(snpMat) <- 1e6+(1:nsnps)
rownames(snpMat) <- 1:n

#generate 3 observations per person, with random intercept:
id <- as.numeric(gl(n,3))
y <- rnorm(3*n) + rnorm(n,sd=.3)[id]
data <- data.frame(y=y, x1=rnorm(3*n,mean=y), x2 =rnorm(3*n), site=sample(6,n,replace=TRUE)[id], id=id )

# set up for GEE with AR-1 correlation structure:
init1 <- boss.set(y~x1+x2+factor(site), id=data$id, data=data, type = "gee", method = "chol", corstr= "ar1")

# set up for linear mixed model (lmm) with random intercept:
init2 <- boss.set(y~x1+x2+factor(site) + (1|id), id=data$id, data=data, type = "lmm", method = "swap")

#fit the model (at one locus):
g <- snpMat[,1]

#compare GEE (ar1) with a linear mixed model (random intercept)
boss.fit(g[init1$id],init1)
#vs
boss.fit(g[init2$id], init2)

#######Using the ncdf wrapper function for all loci:

#####create a net cdf file (see ncdf package)
snpdim <- dim.def.ncdf("position","bases", as.numeric(colnames(snpMat) ))
sampledim <- dim.def.ncdf("sample", "count", as.numeric(rownames(snpMat) ),unlim=TRUE)

varsnp <- var.def.ncdf("snp", "rs", dim = snpdim, missval = -1, prec = "integer")
vargeno <- var.def.ncdf("genotype","base", dim = list(snpdim, sampledim), missval = -1, prec = "double")

ncloc <- "Example.nc"
genofile <- create.ncdf(ncloc,list(varsnp,vargeno))
put.var.ncdf(genofile,"genotype",t(snpMat),start=c(1,1),count=c(-1,-1))
put.var.ncdf(genofile,"genotype",t(snpMat),start=c(1,1),count=c(-1,-1))
#####

nc <- open.ncdf(ncloc)

g.labels <- "position"
id.labels <- "sample"

#gives robust results 
results1 <- boss.ncdf(nc,init1,id.labels,g.labels, gdim = 1 ,chunk=100, robust = TRUE)
results2 <- boss.ncdf(nc,init2,id.labels,g.labels, gdim = 1 ,chunk=100)

}
