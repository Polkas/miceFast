\name{Rcpp_corrData-class}
\Rdversion{1.1}
\docType{class}
\alias{Rcpp_corrData-class}
\alias{corrData}
\title{Class \code{"Rcpp_corrData"}}
\description{
This C++ class could be used to build a corrData object by invoking \code{new(corrData,...)} function.
}

\references{
See the documentation for Rcpp modules for more details of how this class was built.
\code{vignette("Rcpp-modules", package = "Rcpp")}
}

\note{
This is only frame for building C++ object which could be used to implement certain methods.
Check the vignette for more details of implementing methods.

Vigniette:
\url{https://CRAN.R-project.org/package=miceFast}
}

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
