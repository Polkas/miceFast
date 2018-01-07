\name{Rcpp_miceFast-class}
\Rdversion{1.1}
\docType{class}
\alias{Rcpp_miceFast-class}
\alias{miceFast}
\title{Class \code{"Rcpp_miceFast"}}

\description{
This C++ class could be used to build a miceFast objects by invoking \code{new(miceFast,...)} function.
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
#showClass("Rcpp_miceFast")
show(miceFast)
new(miceFast)
}
\keyword{classes}
\section{Methods}{
  \describe{
    \item{\code{initialize(...)}:}{ ~~ }
    \item{\code{finalize()}:}{ ~~ }
    \item{\code{sortData_byg()}:}{ sort data by a grouping variable }
    \item{\code{impute(...)}:}{ impute data }
    \item{\code{get_models(...)}:}{  possible quantitative models for a certain type of dependent variable }
    \item{\code{get_model(...)}:}{  a recommended quantitative model for a certain type of dependent variable }
    \item{\code{get_index_full_R(...)}:}{ get index of rows where dependent and independent variable exists }
    \item{\code{get_index_NA_R(...)}:}{ get index of rows where dependent is missing and independent variable exists }
    \item{\code{is_updated(...)}:}{ is variables at the object modified by impute with a force parameter set to TRUE }
  }
}
