\name{Second_growth_forests}
\alias{Second_growth_forests}
\docType{data}
\title{
  Second-growth forests data (datatype = "incidence_raw" for taxonomic diversity)
}
\description{
  This dataset includes tree incidence data in 100 subplots (each with 0.01 ha) collected from two 
  second-growth forests, namely Cuatro Rios (CR) and Juan Enriquez (JE) in Costa Rica.
  Each 1-ha forest was divided into 100 subplots (each with 0.01 ha) and only species' incidence records 
  in each subplot were used to compute the incidence frequency for a species 
  (i.e., the number of subplots in which that species occurred); see Chazdon et al. (2021, 2022)
  for sampling details and pertinent analyses. The original time series data covers 2005 to 2017; here only two-year data
   (2005 and 2017) are used for assessing temporal beta diversity between the two years within each forest. 
  

}
\usage{data(Second_growth_forests)}
\format{
  Second_growth_forests is a list with two forests. The input format for each forests is a list with two species-by-sampling-units matrices ("Year_2005" and "Year_2017"). Each matrix record the species as 0 (undetect) or 1 (detect) in each sampling units. \cr
  A list of 2  \cr
  $ CR 2005 vs. 2017: A list of 2 \cr
  .. ..$ Year_2005 (151 (species) x 100 (quadrats)) \cr
  .. ..$ Year_2017 (151 (species) x 100 (quadrats)) \cr 
  $ JE 2005 vs. 2017: A list of 2 \cr
  .. ..$ Year_2005 (107 (species) x 100 (quadrats)) \cr
  .. ..$ Year_2017 (107 (species) x 100 (quadrats)) \cr 
}
\references{
  Chazdon, R. (2021). Tree abundance in eight 1-ha tropical forest plots in northeastern Costa Rica from 1997-2017, https://doi.org/10.5061/dryad.ncjsxksvr \cr\cr
  Chazdon, R. L., N. Norden, R. K. Colwell, and A. Chao. (2022). Monitoring recovery of tree diversity during tropical forest restoration: lessons from long-term trajectories of natural regeneration. Philosophical Transactions of the Royal Society B, 378: 20210069. \cr\cr
}

\keyword{datasets}
