\name{planLTPD}
\alias{planLTPD}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculation of LTPD plan for sampling inspection by variables
}
\description{
Calculation of LTPD plan (sample size \eqn{n} and critical value \eqn{k}) for sampling inspection by variables
}
\usage{
planLTPD(N, pt, pbar, b = 0.1, cm = 1, method = c("exact", 
"napprox"), intdif = 20)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{N}{
lot size (number of items in the lot)
}
  \item{pt}{
lot tolerance fraction defective
}
  \item{pbar}{
process average fraction defective
}
  \item{b}{
probability of accepting a lot of tolerance quality \code{pt} (consumer risk), default value 0.1
}
  \item{cm}{
parameter used in cost function of plans (see \code{\link{Ims}}) and defined as cost of inspection of one item by variables divided by cost of inspection of the item by attributes, default value 1
}
  \item{method}{
type of \code{\link{OC}} (operating characteristic) used, either \code{"exact"} (default) using non-central t-distribution, or \code{"napprox"} for approximate OC
}
  \item{intdif}{
parameter used in finding \eqn{n} when \code{method="exact"}, controlling width of interval to be searched in for sampling plan optimal w.r.t. cost, difference from interval centre (\eqn{n} resulting from \code{planLTPD(..., method = "napprox"}) to lower or upper bound, default value 20
}
}
\value{
An instance of \code{\link{ACSPlan-class}}, with sample size in slot \code{n} and critical value in slot \code{k}.
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Klufa, J.: Exact calculation of the Dodge-Romig LTPD single sampling plans for inspection by variables. Statistical Papers, Springer, Vol. 51(2), pages 297-305, 2010.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{LTPDvar-package}}, \code{\link{OC}}, \code{\link{ACSPlan-class}}, \code{\link{Ims}}
}
\examples{
# find LTPD plan
planLTPD(N=1000,pt=0.1,pbar=0.001);
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{misc}
% __ONLY ONE__ keyword per line
