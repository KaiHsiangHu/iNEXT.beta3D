\name{Second_growth_forests}
\alias{Second_growth_forests}
\docType{data}
\title{
Second-growth forests data (datatype = "incidence_raw" for taxonomic diversity)
}
\description{
Second_growth_forests data is a yearly tree species data collected from six 1-ha (50 m x 200 m) second-growth rainforests in Costa Rica. The six forests, with their abbreviated names and ages in 2005 (the selected reference time; see below) are: Finca el Bejuco (FEB, 10 years), Juan Enriquez (JE, 10 years), Lindero Sur (LSUR, 20 years), Tirimbina (TIR, 23 years), Lindero el Peje (LEP, 28 years), and Cuatro Rios (CR, 33 years). The species abundance data for all stems larger or equal than 5 cm in DBH were recorded annually from 2005 to 2017; see Chazdon et al. (2022) for sampling details and summary statistics. All data are available through Dryad (Chazdon, 2021).
}
\usage{data(Second_growth_forests)}
\format{
  Second_growth_forests is a list with six forests. The input format for each forests is a list of several lists incidence raw. In each forest, it has twelve lists: "Year_2005_and_2006", "Year_2005_and_2007", ...,and "Year_2005_and_2017" in it and each lists has two species-by-sampling-units matrices ("Year_2005" and "Year_2006" in "Year_2005_and_2006" as example). Each matrix record the species as 0 (undetect) or 1 (detect) in each sampling units. \cr
 A list of 6  \cr
 $ CR 2005 vs. 2017 :List of 2 \cr
 .. ..$ Year_2005 (168 (species) x 100 (quadrats)) \cr
 .. ..$ Year_2017 (168 (species) x 100 (quadrats)) \cr 
 $ JE 2005 vs. 2017 :List of 2 \cr
 .. ..$ Year_2005 (109 (species) x 100 (quadrats)) \cr
 .. ..$ Year_2017 (109 (species) x 100 (quadrats)) \cr 
 }
\references{
Chazdon, R. 2021. Tree abundance in eight 1-ha tropical forest plots in northeastern Costa Rica from 1997-2017, https://doi.org/10.5061/dryad.ncjsxksvr \cr\cr
Chazdon, R. L., N. Norden, R. K. Colwell, and A. Chao. 2022. Monitoring recovery of tree diversity during tropical forest restoration: lessons from long-term trajectories of natural regeneration. Philosophical Transactions of the Royal Society B, 378: 20210069. \cr\cr
}

\keyword{datasets}
