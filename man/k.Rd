\name{k}
\alias{k}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function for extracting critical value
}
\description{
accessor function for extracting critical value from sampling plan
}
\usage{
k(object)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
sampling plan
}
}
\value{
single value
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
 code{\link{ACSPlan-class}}, 
}
\examples{
# first create an acceptance sampling plan
planek=new("ACSPlan",n=100,k=3)
k(planek)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{misc}
% __ONLY ONE__ keyword per line
