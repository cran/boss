\name{boss.ncdf}
\Rdversion{1.1}
\alias{boss.ncdf}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Wrapper function for anlayzing data in net CDF files
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Takes an object from boss.set() and a net CDF file, reads the data into memory in chunks and performs the analysis
}
\usage{
boss.ncdf(nc, init, id.labels = NULL, g.labels = NULL, subset = NULL, gdim = 1, chunk = 1000, verbose = TRUE, outfile = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nc}{
%%     ~~Describe \code{nc} here~~
a net CDF file
}
  \item{init}{
%%     ~~Describe \code{init} here~~
An object from boss.set()
}
  \item{id.labels}{
%%     ~~Describe \code{id.labels} here~~
Name in the net CDF file where sample id's are given, to be matched with the id's used in boss.set().
}
  \item{g.labels}{
%%     ~~Describe \code{g.labels} here~~
Name in the net CDF file where genotype labels can be found.  Used to label output 
}
  \item{subset}{
%%     ~~Describe \code{subset} here~~
optional, labels of a subset of markers on which to perform the analysis
}
  \item{gdim}{
%%     ~~Describe \code{gdim} here~~
Dimension in the net CDF file where genotype can be found, either 1 (rows) or 2 (columns)
}
  \item{chunk}{
%%     ~~Describe \code{chunk} here~~
Size of chunks to read from the net CDF file.  For 'swap' method, a few thousand at a time is the fastest.  For the 'chol' method, 100 or so gives good efficiency.
}
  \item{verbose}{
%%     ~~Describe \code{verbose} here~~
logical.  Whether or not print progress indicators.  a printed `!' indicates a failed regression, typically due to a constant genotype.
}
  \item{outfile}{
%%     ~~Describe \code{outfile} here~~
optional location of where to write the output, as a csv file
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
aditional arguments to be given to boss.fit()
}
}
\value{
Returns matrix of results giving the MAF of genotype, genotype-related coefficients and variances, and Chi Squared statistics.  If `outfile` is specified, it prints the results to `outfile.csv`, and no object is returned.

}
\references{
\emph{Fast Computation for Genome Wide Association Studies using Boosted One-Step Statistics}
Arend Voorman; Ken Rice; Thomas Lumley
Bioinformatics 2012; doi: 10.1093/bioinformatics/bts291
}

\author{
Arend Voorman
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{boss.fit}}
\code{\link{boss.set}}
\code{\link{ncdf}}
}
\examples{

#generate fake data:
n <- 500
nsnps <- 1000
snpMat <- replicate(nsnps,rbinom(n,2,.2))
colnames(snpMat) <- 1e6+(1:nsnps)
rownames(snpMat) <- 1:n

snpdim <- dim.def.ncdf("position","bases", as.numeric(colnames(snpMat) ))
sampledim <- dim.def.ncdf("sample", "count", as.numeric(rownames(snpMat) ))

varsnp <- var.def.ncdf("snp", "rs", dim = snpdim, missval = -1, prec = "integer")
vargeno <- var.def.ncdf("genotype","base", dim = list(snpdim, sampledim), missval = -1, prec = "double")
genofile <- create.ncdf("Example.nc",list(varsnp,vargeno))

for(i in 1:n){
	put.var.ncdf(genofile,"genotype", snpMat[i,], start=c(1,i), count = c(-1,1))
	}

nc <- open.ncdf("Example.nc")

#generate 3 observations per person, with random intercept:
id <- as.numeric(gl(3*n,3))
y <- rnorm(3*n) + rnorm(n,sd=.3)[id]
data <- data.frame(y=y, x1=rnorm(3*n,mean=y), x2 =rnorm(3*n), site=sample(6,n,replace=TRUE)[id], id=id )

init <- boss.set(y~x1+x2+factor(site), id=data$id, data=data, type = "gee", method = "chol", corstr= "ar1")

g.labels <- "position"
id.labels <- "sample"

results <- boss.ncdf(nc, init, id.labels, g.labels, gdim = 1,chunk=100, robust = TRUE, thresh=0)

close.ncdf(nc)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
