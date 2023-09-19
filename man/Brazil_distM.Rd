\name{Brazil_distM}
\alias{Brazil_distM}
\docType{data}
\title{
Species pairwise distance matrix for \code{Brazil_rainforests} abundance data
}
\description{
This dataset includes a species pairwise symmetric distance matrix for 185 species listed in the 
  \code{Brazil_rainforests} dataset. Each element in the matrix represents a Gower distance between two species
  computed from species traits. \cr
}
\usage{data(Brazil_distM)}
\format{
  Brazil_distM is a 185 x 185 (species by species) symmetric data.frame. The value of Brazil_distM
  is between zero (for identical pairs) and one.  \cr
 $ Carpotroche_brasiliensis      : num  0 0.522 0.522 0.253 0.0.396 ... \cr
 $ Astronium_concinnum           : num  0.522 0 0 0.525 0.0.625 ... \cr
 $ Astronium_graveolens          : num  0.522 0 0 0.525 0.0.625 ... \cr
 ....... \cr
 }


\keyword{datasets}
