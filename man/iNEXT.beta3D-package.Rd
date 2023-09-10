\name{iNEXT.beta3D-package}
\alias{iNEXT.beta3D-package}
\alias{iNEXT.beta3D}
\docType{package}
\title{
iNterpolation and EXTrapolation with beta diversity for three dimensions of biodiversity
}
\description{
iNEXT.beta3D (iNterpolation and EXTrapolation with beta diversity for three dimensions of biodiversity) is a sequel to iNEXT. Here the three dimensions (3D) of diversity means taxonomic diversity (TD), phylogenetic diversity (PD) and functional diversity (FD). A unified framework based on Hill numbers (for TD) and their generalizations (Hill-Chao numbers, for PD and FD) is adopted in the package. TD quantifies the effective number of species, mean-PD (PD divided by tree depth) quantifies the effective number of lineages, and FD quantifies the effective number of virtual functional groups (or functional "species"). Thus, TD, PD (or mean-PD), and FD are all in the same units of species/lineage equivalents and can be meaningfully compared; see Chao et al. (2021) for a review of the unified framework. 

For each of the three dimensions, iNEXT.beta3D focuses on the multiplicative diversity decomposition (alpha, beta and gamma) of orders q = 0, 1 and 2 based on sampling data. Beta diversity quantifies the extent of among-assemblage differentiation, or the changes in species/lineages/functional-groups composition and abundance among assemblages. iNEXT.beta3D features standardized 3D estimates with a common sample size (for alpha and gamma diversity) or sample coverage (for alpha, beta and gamma diversity). iNEXT.beta3D also features standardized estimates of four classes of coverage-based dissimilarity measures. Based on the rarefaction and extrapolation (R/E) method for Hill numbers (TD) for q = 0, 1 and 2, Chao et al. (2023b) developed the pertinent R/E theory for taxonomic beta diversity with applications to real-world spatial, temporal and spatio-temporal data. An application to Gentry's global forest data along with a concise description of the theory is provided in Chao et al. (2023a). The extension to phylogenetic and functional beta diversity is generally parallel.

iNEXT.beta3D features two types of R/E sampling curves:

*	Sample-size-based (or size-based) R/E curves: This type of sampling curve plots standardized 3D gamma and alpha diversity with respect to sample size. Note that the size-based beta diversity is not a statistically valid measure (Chao et al. 2023b) and thus the corresponding sampling curve is not provided. 

*	Coverage-based R/E sampling curves: This type of sampling curve plots standardized 3D gamma, alpha, and beta diversity as well as four classes of dissimilarity measures with respect to sample coverage (an objective measure of sample completeness).

This package contains three main functions: \cr\cr
1. \code{iNEXTbeta3D} computing standardized gamma, alpha, beta diversity, and four dissimilarity-turnover indices for three dimensions: taxonomic, phylogenetic and functional diversity at specified sample coverage or sample size. \cr

2. \code{ggiNEXTbeta3D} visualizing the output from the function \code{iNEXTbeta3D}. \cr

3. \code{DataInfobeta3D} provide basic data information including sample size (or sampling units), observed species richness, estimated sample coverage and other data information for each dimensions. \cr\cr


An Online version of iNEXT.beta3D is also available for users without an R background: \cr
https://chao.shinyapps.io/iNEXT_beta3D/. \cr
In the detailed Online iNEXT.beta3D User's Guide, we illustrate all the running procedures in an easily
accessible way through numerical examples with proper interpretations of portions of the output.
All the data of those illustrative examples are included in this package. \cr\cr
}
\author{
Anne Chao, K.H. Hu \cr

Maintainer: K.H. Hu <river19951224@gmail.com>
}
\references{
*	Chao, A., Chiu, C.-H., Hu, K.-H., and Zeleny, D. (2023a). Revisiting Alwyn H. Gentry's forest transect data: a statistical sampling-model-based approach. Japanese Journal of Statistics and Data Science. https://doi.org/10.1007/s42081-023-00214-1 \cr\cr
*	Chao, A., Henderson, P. A., Chiu, C.-H., Moyes, F., Hu, K.-H., Dornelas, M. and Magurran, A. E. (2021). Measuring temporal change in alpha diversity: a framework integrating taxonomic, phylogenetic and functional diversity and the iNEXT.3D standardization. Methods in Ecology and Evolution, 12, 1926-1940. \cr\cr
*	Chao, A., Thorn, S., Chiu, C.-H., Moyes, F., Hu, K.-H., Chazdon, R. L., Wu, J., Dornelas, M., Zeleny, D., Colwell, R. K., and Magurran, A. E. (2023b). Rarefaction and extrapolation with beta diversity under a framework of Hill numbers: the iNEXT.beta3D standardization. Ecological Monographs, 94, e1588. https://doi.org/10.1002/ecm.1588 \cr\cr
}
\keyword{ package }
