% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{air_miss}
\alias{air_miss}
\title{airquality dataset with additional variables}
\format{
A data frame and data table with 154 observations on 11 variables.
\describe{
\item{Ozone}{	numeric	Ozone (ppb) - Mean ozone in parts per billion from 1300 to 1500 hours at Roosevelt Island}
\item{Solar.R}{	numeric	Solar R (lang) - Solar radiation in Langleys in the frequency band 4000–7700 Angstroms from 0800 to 1200 hours at Central Park}
\item{Wind}{	numeric	Wind (mph) - Average wind speed in miles per hour at 0700 and 1000 hours at LaGuardia Airport}
\item{Temp}{	numeric	Temperature (degrees F) - Maximum daily temperature in degrees Fahrenheit at La Guardia Airport.}
\item{Day}{	numeric	Day of month (1--31)}
\item{Intercept}{ numeric a constant}
\item{index}{ numeric id}
\item{weights}{ numeric positive values weights}
\item{groups}{ factor Month (1--12)}
\item{x_character}{ character discrete version of Solar.R (5-levels)}
\item{Ozone_chac}{ character discrete version of Ozone (7-levels)}
\item{Ozone_f}{ factor discrete version of Ozone (7-levels)}
\item{Ozone_high}{logical Ozone higher than its mean}
}
}
\source{
The data were obtained from the New York State Department of Conservation (ozone data) and the National Weather Service (meteorological data).
}
\usage{
air_miss
}
\description{
airquality dataset with additional variables
}
\details{
Daily readings of the following air quality values for May 1, 1973 (a Tuesday) to September 30, 1973.
}
\examples{
\dontrun{
library(data.table)
data(airquality)
data <- cbind(as.matrix(airquality[, -5]),
  Intercept = 1, index = 1:nrow(airquality),
  # a numeric vector - positive values
  weights = rnorm(nrow(airquality), 1, 0.01),
  # months as groups
  groups = airquality[, 5]
)

# data.table
air_miss <- data.table(data)
air_miss$groups <- factor(air_miss$groups)

# Distribution of Ozone - close to log-normal
# hist(air_miss$Ozone)

# Additional vars
# Make a character variable to show package capabilities
air_miss$x_character <- as.character(cut(air_miss$Solar.R, seq(0, 350, 70)))
# Discrete version of dependent variable
air_miss$Ozone_chac <- as.character(cut(air_miss$Ozone, seq(0, 160, 20)))
air_miss$Ozone_f <- cut(air_miss$Ozone, seq(0, 160, 20))
air_miss$Ozone_high <- air_miss$Ozone > mean(air_miss$Ozone, na.rm = T)
}

}
\references{
Chambers, J. M., Cleveland, W. S., Kleiner, B. and Tukey, P. A. (1983) Graphical Methods for Data Analysis. Belmont, CA: Wadsworth.
}
\keyword{datasets}
