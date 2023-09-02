\name{beetle_inc}
\alias{beetle_inc}
\docType{data}
\title{
Beetles incidence raw data (datatype = "incidence_raw" for TD, PD, and FD)
}
\description{
This data comes from an ecological experiment conducted in the Bavarian Forest National Park in southeastern Germany to study whether felling trees that have fallen due to natural disasters will have an impact on the diversity of the area. The experiment was conducted in a Norway spruce forest covering an area of approximately 1,000 hectares. \cr
The forest was hit by a storm in January 2007, and about 800 hectares of the forest were cut down experimentally and the remaining 200 hectares were left untreated and divided into 22 permanent samples quadrats. Data sources and detailed experimental methods can be found in Thorn et al. (2014) \cr
}
\usage{data(beetle_inc)}
\format{
  beetle_inc is a list with two data.frames in each list. \cr
 $ 2008 vs. 2009 unlogged:List of 2 \cr
  .. .. $ 2008 unlogged: a matrix with 179 species (rows) and 22 sampling units (columns) \cr
  .. .. $ 2009 unlogged: a matrix with 179 species (rows) and 22 sampling units (columns) \cr
 $ 2008 vs. 2009 logged  :List of 2 \cr
  .. .. $ 2008 logged: a matrix with 179 species (rows) and 22 sampling units (columns) \cr
  .. .. $ 2009 logged: a matrix with 179 species (rows) and 22 sampling units (columns) \cr
 }
\references{
Thorn, S., Bassler, C., Gottschalk, T., Hothorn, T., Bussler, H., Raffa, K., & Muller, J. (2014). New insights into the consequences of post-windthrow salvage logging revealed by functional structure of saproxylic beetles assemblages. PloS one, 9(7), e101757.
}

\keyword{datasets}
