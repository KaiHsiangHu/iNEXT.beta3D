\name{Brazil_PDFD_data}
\alias{Brazil_PDFD_data}
\docType{data}
\title{
Brazil data (datatype = "abundance" for PD and FD)
}
\description{
Brazil_PDFD_data is the tree species data collected by Magnago and colleagues between 2011 and 2012 from rain-forest fragments in Espirito Santo State, Brazil (Magnago et al., 2014). Fragmentation of tropical rain forests is one of the major drivers to the loss of global biodiversity. Thus, it is important to understand how species composition changes between fragment edges and interiors. We selected 16 transects (8 Edge transects and the corresponding 8 Interior transects) from fragments with area > 100 ha. \cr
Each transect comprised ten 10 x 10 m plots at 20-m intervals. Each Edge transect was placed about 5 m inside the fragment and parallel to the forest edge, and the corresponding Interior transect was located at least 300 m from the nearest edge. Within each transect, every living tree with a diameter at breast height (DBH) > 4.8 cm and 1.3 m height was recorded and identified to species. The goal here is to assesss how beta diversity between the Edge habitat and the Interior habitat varies with fragment area/size. \cr
Brazil_PDFD_data is used to run examples for taxonomic diversity.
}
\usage{data(Brazil_PDFD_data)}
\format{
  Brazil_PDFD_data is a list with eight species-by-site (two columns(sites): Edge and Interior) data.frames/matrices. \cr
 A list of 8 matrices: \cr
 $ Caliman is a matrix with 148 species(rows) and 2 columns ("Edge" and "Interior"). \cr
 $ Fibria(g) is a matrix with 146 species(rows) and 2 columns ("Edge" and "Interior"). \cr
 $ Marim is a matrix with 119 species(rows) and 2 columns ("Edge" and "Interior"). \cr
 $ Rebio1 is a matrix with 130 species(rows) and 2 columns ("Edge" and "Interior"). \cr
 $ Rebio2 is a matrix with 118 species(rows) and 2 columns ("Edge" and "Interior"). \cr
 $ Rebio3 is a matrix with 111 species(rows) and 2 columns ("Edge" and "Interior"). \cr
 $ RNV(c) is a matrix with 142 species(rows) and 2 columns ("Edge" and "Interior"). \cr
 $ Rochedo is a matrix with 123 species(rows) and 2 columns ("Edge" and "Interior"). \cr
 The species name of Brazil_PDFD_data is matched with tip labels of Brazil_tree and row names (column names) of Brazil_distM.
 }
\references{
Magnago, L. F. S., Edwards, D. P., Edwards, F. A., Magrach, A., Martins, S. V., & Laurance, W. F. 2014. Functional attributes change but functional richness is unchanged after fragmentation of Brazilian Atlantic forests. Journal of Ecology 102: 475-485. \cr
Magnago, L. F. S., Magrach, A., Barlow, J., Schaefer, C. E. G. R., Laurance, W. F., Martins, S. V., & Edwards, D. P. 2017. Do fragment size and edge effects predict carbon stocks in trees and lianas in tropical forests? Functional Ecology 31: 542-552. \cr
Magnago, L. F. S., Rocha, M. F., Meyer, L., Martins, S. v., & Meira-Neto, J. A. A. 2015. Microclimatic conditions at forest edges have significant impacts on vegetation structure in large Atlantic forest fragments. Biodiversity and Conservation 24: 2305-2318. 
}

\keyword{datasets}
