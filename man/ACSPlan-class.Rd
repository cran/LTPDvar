\name{ACSPlan-class}
\Rdversion{1.1}
\docType{class}
\alias{ACSPlan-class}
%\alias{k,ACSPlan-method}
%\alias{n,ACSPlan-method}
%\alias{plot,ACSPlan-method}

\title{Class \code{ACSPlan}}
\description{
Class for single-sample plan of sampling inspection by variables. The plan is specified by sample size \eqn{n} and critical value \eqn{k}.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("ACSPlan", ...)}.
Objects represent sampling plan.
}
\section{Slots}{
  \describe{
    \item{\code{n}:}{Object of class \code{"numeric"}, sample size, i. e. number of items to be inspected }
    \item{\code{k}:}{Object of class \code{"numeric"}, critical value }
  }
}
\section{Methods}{
  \describe{
    \item{k}{\code{signature(object = "ACSPlan")}: accessor function for extraction of critical value of the sampling plan }
    \item{n}{\code{signature(object = "ACSPlan")}: accessor function for extraction of sample size of the sampling plan }
    \item{plot}{\code{signature(x = "ACSPlan")}: function for getting operating characteristics plot of the sampling plan }
	 }
}
\references{
%%  ~~put references to the literature/web site here~~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%%  ~~objects to See Also as \code{\link{~~fun~~}}, ~~~
%%  ~~or \code{\linkS4class{CLASSNAME}} for links to other classes ~~~
}
\examples{
showClass("ACSPlan")
}
\keyword{classes}
