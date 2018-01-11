\name{Rcpp_miceFast-class}
\Rdversion{1.1}
\docType{class}
\alias{Rcpp_miceFast-class}
\alias{miceFast}
\title{Class \code{"Rcpp_miceFast"}}

\description{
This C++ class could be used to build a miceFast objects by invoking \code{new(miceFast)} function.
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
Check the vignette for more details of implementing these methods.

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
    \item{\code{set_data(...)}:}{ providing the data by a reference }
    \item{\code{set_g(...)}:}{ providing a grouping variable }
    \item{\code{set_w(...)}:}{ providing a weightinh variable by a reference }
    \item{\code{impute(...)}:}{ impute data under characterstics from the object like a optional grouping or weighting variable }
    \item{\code{impute_force(...)}:}{ impute data under characterstics from the object like a optional grouping or weighting variable and permanently update the variable at the object. WARNING, use it only if you are sure about model parameters. }
    \item{\code{get_models(...)}:}{  possible quantitative models for a certain type of dependent variable }
    \item{\code{get_model(...)}:}{  a recommended quantitative model for a certain type of dependent variable }
    \item{\code{which_updated(...)}:}{ which variables at the object was modified by impute_force}
    \item{\code{initialize(...)}:}{ ... }
    \item{\code{finalize()}:}{ ... }
  }
}
