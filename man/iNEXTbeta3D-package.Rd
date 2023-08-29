\name{iNEXT.beta3D-package}
\alias{iNEXT.beta3D-package}
\alias{iNEXT.beta3D}
\docType{package}
\title{
iNterpolation and EXTrapolation with beta diversity for three dimensions of diversity
}
\description{
iNEXT.beta3D (iNterpolation and EXTrapolation with beta diversity for three dimensions of diversity) is a sequel to iNEXT Online. Here three dimensions (3D) of diversity means taxonomic diversity (TD), phylogenetic diversity (PD) and functional diversity (FD). For each of the three dimensions, iNEXT.beta3D features standardized 3D estimates with a common sample size (for alpha and gamma diversity) or sample coverage (for alpha, beta and gamma diversity). iNEXT.beta3D also features standardized estimates of four classes of coverage-based dissimilarity measures. Based on the rarefaction and extrapolation (R/E) method for Hill numbers (taxonomic diversity) for q = 0, 1 and 2, Chao et al. (2023b) developed the pertinent R/E theory for taxonomic beta diversity with applications to real-world spatial, temporal and spatio-temporal data. An application to Gentry's global forest data along with a concise description of the theory is provided in Chao et al. (2023a). The extension to PD and FD (both are based on Hill-Chao numbers) is generally parallel. 
iNEXT.beta3D features two types of R/E sampling curves:

*	Sample-size-based (or size-based) R/E curves: This type of sampling curve plots standardized 3D gamma and alpha diversity with respect to sample size. Note that the size-based beta diversity is not a statistically valid measure (Chao et al. 2023b) and thus the corresponding sampling curve is not provided. 

*	Coverage-based R/E sampling curves: This type of sampling curve plots standardized 3D gamma, alpha, and beta diversity as well as four classes of dissimilarity measures with respect to sample coverage (an objective measure of sample completeness).

This package contains three main functions: \cr\cr
1. \code{iNEXTbeta3D} Computing standardized gamma, alpha, beta diversity, and four dissimilarity-turnover indices for three dimensions: taxonomic, phylogenetic and functional diversity at specified sample coverage or sample size. \cr\cr
2. \code{ggiNEXTbeta3D} Visualizing the output from the function `iNEXTbeta3D`. \cr\cr
3. \code{DataInfobeta} Provide basic data information including sample size (or sampling units), observed species richness, estimated sample coverage and other data information for each dimensions. \cr\cr

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
Chao, A., Chiu, C.-H., Hu, K.-H., and Zeleny, D. (2023a). Revisiting Alwyn H. Gentry's forest transect data: a statistical sampling-model-based approach. Japanese Journal of Statistics and Data Science. https://doi.org/10.1007/s42081-023-00214-1 \cr\cr
Chao, A., Thorn, S., Chiu, C.-H., Moyes, F., Hu, K.-H., Chazdon, R. L., Wu, J., Dornelas, M., Zeleny, D., Colwell, R. K., and Magurran, A. E. (2023b). Rarefaction and extrapolation with beta diversity under a framework of Hill numbers: the iNEXT.beta3D standardization. Ecological Monographs, 94, e1588. https://doi.org/10.1002/ecm.1588 \cr\cr
Chao, A. and Hu, K.-H. (2023). iNEXT.beta3D Online: interpolation and extrapolation with beta diversity for three dimensions of diversity. Code and user's guide available from http://chao.stat.nthu.edu.tw/wordpress/software_download/inext-beta3d/  \cr\cr
}
\keyword{ package }
