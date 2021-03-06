\name{planAOQL}
\alias{planAOQL}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculation of AOQL plan for sampling inspection by variables
}
\description{
Calculation of AOQL plan (sample size \eqn{n} and critical value \eqn{k}) for sampling inspection by variables. Plans minimize mean inspection cost per lot of process average quality and at the same time satisfy limit on average outgoing quality (see \code{\link{AOQ}}).
}
\usage{
planAOQL(N, pbar, pL, method = c("exact", "napprox","ewmaSK","ewma2"), cm = 1,
	intdif = 20,lam=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{N}{
lot size (number of items in the lot)
}
  \item{pbar}{
process average fraction defective
}
  \item{pL}{
average outgoing quality limit 
}
  \item{method}{
type of \code{\link{OC}} (operating characteristic) used, either \code{"exact"} (default) using non-central t-distribution, or \code{"napprox"} for approximate OC; if EWMA statistic is to be used, \code{"ewmaSK"} stands for the sigma known case,  \code{"ewma2"} stands for the sigma unknown case
}
  \item{cm}{
parameter used in cost function of plans (see \code{\link{Ims}}) and defined as cost of inspection of one item by variables divided by cost of inspection of the item by attributes, default value 1
}
  \item{intdif}{
parameter used in finding \eqn{n} when \code{method} other than \code{"napprox"}, controlling the width of interval to be searched in for sampling plan optimal w.r.t. cost, difference from interval centre (\eqn{n} resulting from \code{planLTPD(..., method = "napprox"}) to lower or upper bound, default value 20
}
 \item{lam}{
smoothing parameter in case that EWMA statistic is used
}
}
\value{
\code{\link{ACSPlan-class}} object
}
\references{
Klufa J (1997) Dodge-Romig AOQL single sampling plans for inspection
by variables. Statistical Papers 38: 111 - 119
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{LTPDvar-package}}, \code{\link{OC}}, \code{\link{AOQ}}, \code{\link{ACSPlan-class}}, \code{\link{Ims}}
}
\examples{
# find AOQL plan
planAOQL(N=1000,pbar=0.005,pL=0.01, method="napprox", cm=1.5);
planAOQL(N=8000, pbar=0.003, pL=0.01, cm=1.5,method="ewmaSK", lam=0.9,intdif=40);
planAOQL(N=8000, pbar=0.003, pL=0.01, cm=1.5,method="ewma2", lam=0.9);
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{misc}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
