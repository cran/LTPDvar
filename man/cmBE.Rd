\name{cmBE}
\alias{cmBE}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Break-even value of cm parameter
}
\description{Break-even value of \eqn{c_m} parameter (which is ratio of cost
of inspection of one item by variables to cost of inspection of the item by
attributes), i. e. the value of  \eqn{c_m} for which mean inspection cost per lot of process average quality for inspection by variables and attributes is
equal to mean inspection cost per lot of process average quality for inspection
by attributes, using plan \eqn{(n, c)}.
}
\usage{
cmBE(N, pbar,px,n,c,type=c("LTPD","AOQL"),
type2 = c("exact", "napprox","ewmaSK","ewma2"),lam=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{N}{
lot size (number of items in the lot)
}
  \item{pbar}{
process average fraction defective
}
  \item{px}{
lot tolerance fraction defective \eqn{pt} (if type="LTPD", which is default) or average outgoing quality limit \eqn{p_L} (if type="AOQL") 
}
  \item{n}{
sample size of benchmark plan for sampling inspection by attributes
}
  \item{c}{
acceptance number of benchmark plan for sampling inspection by attributes
}
  \item{type}{
type of acceptance sampling plan; \code{"LTPD"} stands for Lot Tolerance Percent Defective (see \code{\link{planLTPD}}), \code{"AOQL"} stands for  Average Outgoing Quality Limit (see \code{\link{planAOQL}})
}
  \item{type2}{
type of  {\link{OC}} to be used
}
 \item{lam}{
smoothing parameter in case that the EWMA statistic is to be used, defaults to 1
}
}
\value{
single number
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Kasprikova, N. and Klufa, J.: AOQL Sampling Plans for Inspection by Variables and Attributes Versus the Plans for Inspection by Attributes. Quality Technology & Quantitative Management, 12/6. 2015.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{planLTPD}},  \code{\link{planAOQL}}
}
\examples{
cmBE(N=1000,pbar=0.001,px=0.01,n=80,c=0,type="LTPD",type2="exact");
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{misc}
% __ONLY ONE__ keyword per line
