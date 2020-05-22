\name{Rcpp_corrData-class}
\Rdversion{1.1}
\docType{class}
\alias{Rcpp_corrData-class}
\alias{corrData}
\title{Class \code{"Rcpp_corrData"}}
\description{
This C++ class could be used to build a corrData object by invoking \code{new(corrData,...)} function.
}
\section{Extends}{
Class \code{"\linkS4class{C++Object}"}, directly.

All reference classes extend and inherit methods from \code{"\linkS4class{envRefClass}"}.

}
\references{
See the documentation for RcppArmadillo and Rcpp for more details of how this class was built.
}

\note{
This is only frame for building C++ object which could be used to implement certain methods.
Check the vignette for more details of implementing methods.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
#showClass("Rcpp_corrData")
show(corrData)
}
\keyword{classes}
\section{Methods}{
  \describe{
    \item{\code{initialize(...)}:}{ ~~ }
    \item{\code{finalize()}:}{ ~~ }
    \item{\code{fill(...)}:}{ generating data }
  }
}
