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
Ims(n, N, pt, pbar, cm = 1, b = 0.1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
sample size
}
  \item{N}{
lot size (number of items in the lot)
}
  \item{pt}{
lot tolerance fraction defective
}
  \item{pbar}{
process average fraction defective
}
  \item{cm}{
cost of inspection of one item by variables divided by cost of inspection of the item by attributes, default value 1
}
  \item{b}{
probability of accepting a lot of tolerance quality \code{pt} (consumer risk), default value 0.1
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
Ims(16,100,0.01,0.001)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{misc}
% __ONLY ONE__ keyword per line
