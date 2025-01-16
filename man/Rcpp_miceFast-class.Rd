\name{Rcpp_miceFast-class}
\Rdversion{1.1}
\docType{class}
\alias{Rcpp_miceFast-class}
\alias{miceFast}
\title{Class \code{"Rcpp_miceFast"}}

\description{
This C++ class could be used to build a miceFast objects by invoking \code{new(miceFast)} function.
}

\references{
See the documentation for Rcpp modules for more details of how this class was built.
\code{vignette("Rcpp-modules", package = "Rcpp")}
}

\note{
This is only frame for building C++ object which could be used to implement certain methods.
Check the vignette for more details of implementing these methods.

Vigniette:
\url{https://CRAN.R-project.org/package=miceFast}
}


\examples{
#showClass("Rcpp_miceFast")
show(miceFast)
new(miceFast)
}
\keyword{classes}
\section{Methods}{
  \describe{
    \item{\code{set_data(...)}:}{ providing data by a reference - a numeric matrix}
    \item{\code{set_g(...)}:}{ providing a grouping variable by a reference - a numeric vector WITOUT NA values - positive values}
    \item{\code{set_w(...)}:}{ providing a weightinh variable by a reference - a numeric vector WITOUT NA values - positive values}
    \item{\code{set_ridge(...)}:}{ providing a ridge i.e. the disturbance to diag of XX, default 1e-6}
    \item{\code{get_data(...)}:}{ retrieving the data }
    \item{\code{get_w(...)}:}{ retrieving the weighting variable }
    \item{\code{get_g(...)}:}{ retireiving the grouping variable }
    \item{\code{get_ridge(...)}:}{ retireiving the ridge disturbance }
    \item{\code{get_index(...)}:}{ getting the index }
    \item{\code{impute(...)}:}{ impute data under characterstics from the object like a optional grouping or weighting variable }
    \item{\code{impute_N(...)}:}{ multiple imputations - impute data under characterstics from the object like a optional grouping or weighting variable }
    \item{\code{update_var(...)}:}{ permanently update the variable at the object and data. Use it only if you are sure about model parameters}
    \item{\code{get_models(...)}:}{  get possible quantitative models for a certain type of dependent variable }
    \item{\code{get_model(...)}:}{  get a recommended quantitative model for a certain type of dependent variable }
    \item{\code{which_updated(...)}:}{ which variables at the object was modified by update_var}
    \item{\code{sort_byg(...)}:}{ sort data by the grouping variable }
    \item{\code{is_sorted_byg(...)}:}{check if data is sorted by the grouping variable }
    \item{\code{vifs(...)}:}{Variance inflation factors (VIF) - helps to check when the predictor variables are not linearly related}
    \item{\code{initialize(...)}:}{ ... }
    \item{\code{finalize()}:}{ ... }
  }
}
