\name{Ims}
\alias{Ims}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Inspection cost function
}
\description{
mean inspection cost per lot of process average quality, assuming that the sample is inspected by variables and the remainder of rejected lots is inspected by attributes, divided by parameter \code{cm} (cost of inspecting one item by variables divided by cost of inspecting the item by attributes)
}
\usage{
Ims(n, k, N,  pbar, cm = 1,type = c("exact", "napprox","ewmaSK","ewma2"),lam=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
sample size
}
  \item{k}{
critical value of the samping plan
}
  \item{N}{
lot size (number of items in the lot)
}
  \item{pbar}{
process average fraction defective
}
  \item{cm}{
cost of inspection of one item by variables divided by cost of inspection of the item by attributes, default value 1
}
 \item{type}{
type of \code{\link{OC}} (operating characteristic) used, either \code{"exact"} (default) using non-central t-distribution, or \code{"napprox"} for approximate OC; if EWMA statistic is to be used, \code{"ewmaSK"} stands for the sigma known case,  \code{"ewma2"} stands for the sigma unknown case
}
 \item{lam}{
smoothing parameter in case that EWMA statistic is used
}
}
\value{
single value 
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{LTPDvar-package}}
}
\examples{
Ims(20, 2.58555,1000, 0.001 ,1.5,type="ewmaSK",lam=1 )
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{misc}
% __ONLY ONE__ keyword per line
