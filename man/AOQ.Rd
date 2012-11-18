\name{AOQ}
\alias{AOQ}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Average outgoing quality
}
\description{
Average outgoing quality is the mean fraction defective after inspection when the fraction defective before inspection was \eqn{p}, lot size is \eqn{N} and plan \eqn{(n,k)} is used for sampling inspection. The average outgoing quality (assuming that all defective items found are replaced by good ones) is approximately
\deqn{AOQ(p)=\left(1-\frac{n}N\right)\!\cdot p\cdot\!L(p;n,c).}
}
\usage{
AOQ(p,n,k,N, type=c("t","napprox"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p}{
fraction defective before inspection
}
  \item{n}{
sample size
}
  \item{k}{
critical value
}
  \item{N}{
lot size (number of items in the lot)
}
  \item{type}{
type of operating characteristic, see \code{\link{OC}}, default value "t" i. e. exact OC is used
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
AOQ(0.002,41,2.057083,1000)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{misc}
% __ONLY ONE__ keyword per line
