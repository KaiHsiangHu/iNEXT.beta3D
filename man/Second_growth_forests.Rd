\name{Second_growth_forests}
\alias{Second_growth_forests}
\docType{data}
\title{
Second-growth forests data (datatype = "incidence_raw" for taxonomic diversity)
}
\description{
Second_growth_forests data is a yearly tree species data collected from six 1-ha (50 m x 200 m) second-growth rainforests in Costa Rica. The two forests, with their abbreviated names and ages in 2005 are: Juan Enriquez (JE, 10 years) and Cuatro Rios (CR, 33 years). The species abundance data for all stems larger or equal than 5 cm in DBH were recorded annually from 2005 to 2017; see Chazdon et al. (2022) for sampling details and summary statistics. Here we use year 2005 and 2017 in JE and CR forests to compute iNEXT beta. All data are available through Dryad (Chazdon, 2021).
}
\usage{data(Second_growth_forests)}
\format{
  Second_growth_forests is a list with two forests. The input format for each forests is a list with two species-by-sampling-units matrices ("Year_2005" and "Year_2017"). Each matrix record the species as 0 (undetect) or 1 (detect) in each sampling units. \cr
 A list of 2  \cr
 $ CR 2005 vs. 2017 :A list of 2 \cr
 .. ..$ Year_2005 (168 (species) x 100 (quadrats)) \cr
 .. ..$ Year_2017 (168 (species) x 100 (quadrats)) \cr 
 $ JE 2005 vs. 2017 :A list of 2 \cr
 .. ..$ Year_2005 (109 (species) x 100 (quadrats)) \cr
 .. ..$ Year_2017 (109 (species) x 100 (quadrats)) \cr 
 }
\references{
Chazdon, R. 2021. Tree abundance in eight 1-ha tropical forest plots in northeastern Costa Rica from 1997-2017, https://doi.org/10.5061/dryad.ncjsxksvr \cr\cr
Chazdon, R. L., N. Norden, R. K. Colwell, and A. Chao. 2022. Monitoring recovery of tree diversity during tropical forest restoration: lessons from long-term trajectories of natural regeneration. Philosophical Transactions of the Royal Society B, 378: 20210069. \cr\cr
}

\keyword{datasets}
