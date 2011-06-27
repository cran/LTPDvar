\name{OC}
\alias{OC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Operating characteristics
}
\description{Calculation of probability of acceptance of a lot with fraction defective \eqn{p})
when using plan \eqn{(n, k)} for sampling inspection
}
\usage{
OC(p, n, k, type = c("t", "napprox"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p}{
fraction defective
}
  \item{n}{
sample size
}
  \item{k}{
critical value
}
  \item{type}{
function used for operating characteristics; \code{"t"} stands for exact function using non-central t-distribution (see (Johnson  and Welch, 1940)), \code{"napprox"} stands for approximation using normal distribution,  (see (Jennett  and Welch, 1939))
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
probability of acceptance of a lot (single number)
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Jennett, W. J. -  Welch, B. L.: The Control of Proportion Defective as Judged by a Single Quality Characteristic Varying on a Continuous Scale, Supplement to the Journal of the Royal Statistical Society, Vol. 6, No. 1, pp. 80-88, 1939.

Johnson, N. L. - Welch, B. L.: Applications of the Non-Central t-Distribution, Biometrika, Vol. 31, No. 3/4, pp. 362-389, 1940.
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
OC(p=0.1,n=85,k=2.44)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{htest}
% __ONLY ONE__ keyword per line
