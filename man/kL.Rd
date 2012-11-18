\name{kL}
\alias{kL}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function that finds critical value
}
\description{
Function that finds critical value for LTPD acceptance sampling by variables.
}
\usage{
kL(n, pt, b = 0.1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
sample size
}
  \item{pt}{
lot tolerance fraction defective
}
  \item{b}{
probability of accepting a lot of tolerance quality \code{pt} (consumer risk), default value 0.1
}
}
\value{
critical value (single number)
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{planLTPD}},  
}
\examples{
kL(16,pt=0.1)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{htest}
% __ONLY ONE__ keyword per line
