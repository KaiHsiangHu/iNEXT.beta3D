#' iNterpolation and EXTrapolation with beta diversity for TD, PD and FD 
#' 
#' \code{iNEXTbeta3D}: compute standardized 3D estimates with a common sample size
#'(for alpha and gamma diversity) or sample coverage (for alpha, beta, gamma diversity as well as
#' dissimilarity indices; see Chao et al. (2023) for the theory.
#'                      
#' @param data (a) For \code{datatype = "abundance"}, species abundance data for a single dataset can be input as a \code{matrix/data.frame} (species-by-assemblage); data for multiple datasets can be input as a \code{list} of \code{matrices/data.frames}, with each matrix representing a species-by-assemblage abundance matrix for one of the datasets.\cr
#' (b) For \code{datatype = "incidence_raw"}, data for a single dataset with N assemblages can be input as a \code{list} of \code{matrices/data.frames}, with each matrix representing a species-by-sampling-unit incidence matrix for one of the assemblages; data for multiple datasets can be input as multiple lists.
#' @param diversity selection of diversity type: \code{'TD'} = Taxonomic diversity, \code{'PD'} = Phylogenetic diversity, and \code{'FD'} = Functional diversity.
#' @param q a numerical vector specifying the diversity orders. Default is \code{q = c(0, 1, 2)}.
#' @param datatype data type of input data: individual-based abundance data (\code{datatype = "abundance"}) or species by sampling-units incidence matrix (\code{datatype = "incidence_raw"}) with all entries being \code{0} (non-detection) or \code{1} (detection).
#' @param base standardization base: coverage-based rarefaction and extrapolation for gamma, alpha, beta diversity, and four classes of dissimilarity indices (\code{base = "coverage"}), or sized-based rarefaction and extrapolation for gamma and alpha diversity (\code{base = "size"}). Default is \code{base = "coverage"}.
#' @param level a numerical vector specifying the particular values of sample coverage (between 0 and 1 when 
#' \code{base = "coverage"}) or sample sizes (\code{base = "size"}) that will be used to compute standardized
#'  diversity/dissimilarity. Asymptotic diversity estimator can be obtained by setting \code{level = 1}
#'   ( i.e., complete coverage for \code{base = "coverage"}). By default (with \code{base = "coverage"}), 
#'   this function computes the standardized 3D gamma, alpha, beta diversity, and four dissimilarity indices 
#'   for coverage up to one (for \code{q = 1, 2}) or up to the coverage of double the reference sample size 
#'   (for \code{q = 0}), in increments of 0.025. The extrapolation limit for beta diversity is defined as that
#'    for alpha diversity. If users set \code{base = "size"}, this function computes the size-based standardized
#'     3D gamma and alpha diversity estimates based on 40 equally-spaced sample sizes/knots from sample size 1 
#'     up to double the reference sample size.
#' @param nboot a positive integer specifying the number of bootstrap replications when assessing sampling uncertainty and constructing confidence intervals. Bootstrap replications are generally time consuming. Enter \code{0} to skip the bootstrap procedures. Default is \code{nboot = 10}. If more accurate results are required, set \code{nboot = 100} (or \code{200}).
#' @param conf a positive number < 1 specifying the level of confidence interval. Default is \code{conf = 0.95}.
#' @param PDtree (required argument for \code{diversity = "PD"}), a phylogenetic tree in Newick format for all observed species in the pooled dataset. 
#' @param PDreftime (argument for \code{diversity = "PD"}), a numerical value specifying reference time for PD. Default is \code{PDreftime = NULL} (i.e., the age of the root of PDtree).  
#' @param PDtype (argument for \code{diversity = "PD"}), select PD type: \code{PDtype = "PD"} (effective total branch length) or \code{PDtype = "meanPD"} (effective number of equally divergent lineages). Default is \code{PDtype = "meanPD"}, where \code{meanPD = PD/tree depth}.
#' @param FDdistM (required argument for \code{diversity = "FD"}), a species pairwise distance matrix for all species in the pooled dataset. 
#' @param FDtype (argument for \code{diversity = "FD"}), select FD type: \code{FDtype = "tau_value"} for FD under a specified threshold value, or \code{FDtype = "AUC"} (area under the curve of tau-profile) for an overall FD which integrates all threshold values between zero and one. Default is \code{FDtype = "AUC"}.  
#' @param FDtau (argument for \code{diversity = "FD"} and \code{FDtype = "tau_value"}), a numerical value between 0 and
#'  1 specifying the tau value (threshold level) that will be used to compute FD. If \code{FDtype = NULL} (default), 
#'  then threshold level is set to be the mean distance between any two individuals randomly selected from the pooled 
#'  dataset (i.e., quadratic entropy). 
#' @param FDcut_number (argument for \code{diversity = "FD"} and \code{FDtype = "AUC"}), a numeric number to cut [0, 1] interval into equal-spaced sub-intervals to obtain the AUC value by integrating the tau-profile. Equivalently, the number of tau values that will be considered to compute the integrated AUC value. Default is \code{FDcut_number = 30}. A larger value can be set to obtain more accurate AUC value.
#' 
#' @import tidyverse
#' @import magrittr
#' @import ggplot2
#' @import abind
#' @import ape
#' @import phytools
#' @import phyclust
#' @import tidytree
#' @import colorRamps
#' @import iNEXT.3D
#' @import future.apply
#' @import ade4
#' @import tidyr
#' @import tibble
#' 
#' @return For \code{base = "coverage"}, return a list of seven data frames with three diversity (gamma, alpha, and beta
#'  diversity) and four dissimilarity measures. For \code{base = "size"}, return a list of two matrices with two diversity
#'  (gamma and alpha diversity).\cr 
#'  
#'  For \code{base = "coverage"}, the output in each data frame includes: 
#'  \item{Dataset}{the name of dataset. }
#'  \item{Order.q}{the diversity order of q.} 
#'  \item{SC}{the target standardized coverage value.}
#'  \item{Size}{the corresponding sample size.}
#'  \item{Alpha/Beta/Gamma/Dissimilarity}{the estimated diversity/dissimilarity estimate.}
#'  \item{Method}{Rarefaction, Observed, or Extrapolation, depending on whether the target coverage is less than, equal to, or greater than the coverage of the reference sample.} 
#'  \item{s.e.}{standard error of standardized estimate.}
#'  \item{LCL, UCL}{the bootstrap lower and upper confidence limits for the diversity/dissimilarity with a default significance level of 0.95.}
#'  Similar output is obtained for \code{base = "size"}.\cr
#'  
#'  
#' @examples
#' ## Taxonomic diversity for abundance data
#' # Coverage-based standardized TD estimates and related statistics
#' data(Brazil_rainforests)
#' output1c = iNEXTbeta3D(data = Brazil_rainforests, diversity = 'TD', 
#'                        datatype = 'abundance', base = "coverage", nboot = 30, conf = 0.95)
#' output1c
#' 
#' 
#' # Size-based standardized TD estimates and related statistics
#' data(Brazil_rainforests)
#' output1s = iNEXTbeta3D(data = Brazil_rainforests, diversity = 'TD', 
#'                        datatype = 'abundance', base = "size", nboot = 30, conf = 0.95)
#' output1s
#' 
#' 
#' ## Phylogenetic diversity for abundance data
#' # Coverage-based standardized PD estimates and related statistics
#' data(Brazil_rainforests)
#' data(Brazil_tree)
#' output2c = iNEXTbeta3D(data = Brazil_rainforests, diversity = 'PD', 
#'                        datatype = 'abundance', base = "coverage", nboot = 10, conf = 0.95, 
#'                        PDtree = Brazil_tree, PDreftime = NULL, PDtype = 'PD')
#' output2c
#' 
#' # Size-based standardized PD estimates and related statistics
#' data(Brazil_rainforests)
#' data(Brazil_tree)
#' output2s = iNEXTbeta3D(data = Brazil_rainforests, diversity = 'PD', 
#'                        datatype = 'abundance', base = "size", nboot = 20, conf = 0.95, 
#'                        PDtree = Brazil_tree, PDreftime = NULL, PDtype = 'PD')
#' output2s
#' 
#' 
#' ## Functional diversity for abundance data under a specified threshold level
#' # Coverage-based standardized FD estimates and related statistics
#' data(Brazil_rainforests)
#' data(Brazil_distM)
#' output3c = iNEXTbeta3D(data = Brazil_rainforests, diversity = 'FD', 
#'                        datatype = 'abundance', base = "coverage", nboot = 30, conf = 0.95, 
#'                        FDdistM = Brazil_distM, FDtype = 'tau_value', FDtau = NULL)
#' output3c
#' 
#' # Size-based standardized FD estimates and related statistics
#' data(Brazil_rainforests)
#' data(Brazil_distM)
#' output3s = iNEXTbeta3D(data = Brazil_rainforests, diversity = 'FD', 
#'                        datatype = 'abundance', base = "size", nboot = 30, conf = 0.95, 
#'                        FDdistM = Brazil_distM, FDtype = 'tau_value', FDtau = NULL)
#' output3s
#' 
#' 
#' ## Functional diversity for abundance data when all thresholds from 0 to 1 are considered
#' # Coverage-based standardized TD estimates and related statistics
#' data(Brazil_rainforests)
#' data(Brazil_distM)
#' output4c = iNEXTbeta3D(data = Brazil_rainforests, diversity = 'FD', 
#'                        datatype = 'abundance', base = "coverage", nboot = 0, conf = 0.95, 
#'                        FDdistM = Brazil_distM, FDtype = 'AUC', FDcut_number = 30)
#' output4c
#' 
#' 
#' # Size-based standardized TD estimates and related statistics
#' data(Brazil_rainforests)
#' data(Brazil_distM)
#' output4s = iNEXTbeta3D(data = Brazil_rainforests, diversity = 'FD', 
#'                        datatype = 'abundance', base = "size", nboot = 10, conf = 0.95, 
#'                        FDdistM = Brazil_distM, FDtype = 'AUC', FDcut_number = 30)
#' output4s
#' 
#' 
#' ## Taxonomic diversity for incidence data
#' # Coverage-based standardized TD estimates and related statistics
#' data(Second_growth_forests)
#' output5c = iNEXTbeta3D(data = Second_growth_forests, diversity = 'TD', datatype = 'incidence_raw', 
#'                        base = "coverage", level = NULL, nboot = 10, conf = 0.95)
#' output5c
#' 
#' 
#' # Size-based standardized TD estimates and related statistics
#' data(Second_growth_forests)
#' output5s = iNEXTbeta3D(data = Second_growth_forests, diversity = 'TD', datatype = 'incidence_raw', 
#'                        base = "size", level = NULL, nboot = 10, conf = 0.95)
#' output5s
#' 
#' 
#' @references
#' Chao, A., Thorn, S., Chiu, C.-H., Moyes, F., Hu, K.-H., Chazdon, R. L., Wu, J., Dornelas, M., Zeleny, D., Colwell, 
#' R. K., and Magurran, A. E. (2023). Rarefaction and extrapolation with beta diversity under a framework of Hill numbers:
#'  the iNEXT.beta3D standardization. Ecological Monographs e1588.
#' @export
iNEXTbeta3D = function(data, diversity = 'TD', q = c(0, 1, 2), datatype = 'abundance', 
                       base = 'coverage', level = NULL, nboot = 10, conf = 0.95, 
                       PDtree = NULL, PDreftime = NULL, PDtype = 'meanPD',
                       FDdistM = NULL, FDtype = 'AUC', FDtau = NULL, FDcut_number = 30) {
  
  max_alpha_coverage = F
  
  
  ## Check parameter setting
  if (is.na(pmatch(diversity, c("TD", "PD", "FD")))) stop("invalid diversity")
  
  if (!inherits(q, "numeric"))
    stop("invlid class of order q, q should be a postive value/vector of numeric object", call. = FALSE)
  if (min(q) < 0){
    warning("ambigous of order q, we only compute postive q.", call. = FALSE)
    q <- q[q >= 0]
  }
  
  if (datatype == "incidence" | datatype == "incidence_freq") stop('Please try datatype = "incidence_raw".')  
  if(is.na(pmatch(datatype, c("abundance", "incidence_raw"))))
    stop("invalid datatype")
  
  if (is.na(pmatch(base, c("size", "coverage")))) stop("invalid datatype")
  
  if (! (is.null(level) | inherits(level, "numeric")))
    stop("invlid class of level, level should be a postive value/vector of numeric object", call. = FALSE)
  
  
  if ((nboot < 0) | (is.numeric(nboot) == F)) stop('Please enter non-negative integer for nboot.', call. = FALSE)
  
  if ((conf < 0) | (conf > 1) | (is.numeric(conf) == F)) stop('Please enter value between zero and one for confident interval.', call. = FALSE)
  
  if (! (is.null(PDreftime) | inherits(PDreftime, "numeric")))
    stop("invalid class of reference time, PDreftime should be a postive value of numeric object.", call. = FALSE)
  if (length(PDreftime) > 1)
    stop("PDreftime can only accept a value instead of a vector.", call. = FALSE)
  
  if(is.na(pmatch(PDtype, c("PD", "meanPD"))))
    stop("Incorrect type of phylogenetic diversity type, please use either 'PD' or 'meanPD'.", call. = FALSE)
  
  if (FDtype == "tau_values") stop('Please try FDtype = "tau_value".')  
  if(is.na(pmatch(FDtype, c("AUC", "tau_value"))))
    stop("Incorrect type of functional diversity type, please use either 'AUC' or 'tau_value'", call. = FALSE)
  
  if (! (is.null(FDtau) | inherits(FDtau, "numeric")))
    stop("invalid class of tau value, FDtau should be a postive value between zero and one.", call. = FALSE)
  if (length(FDtau) > 1)
    stop("FDtau only accept a value instead of a vector.", call. = FALSE)
  
  if (!inherits(FDcut_number, "numeric"))
    stop("invalid class of FD cut number, FDcut_number should be a postive value.", call. = FALSE)
  if (FDcut_number < 2)
    stop("invalid FDcut_number, FDcut_number should be a postive value larger than one.", call. = FALSE)
  if (length(FDcut_number) > 1)
    stop("FDcut_number only accept a value instead of a vector.", call. = FALSE)
  
  ##
  
  
  
  if (datatype == 'abundance') {
    
    if ( inherits(data, "data.frame") | inherits(data, "matrix") ) data = list(Dataset_1 = data)
    
    if (class(data) == "list"){
      
      if (is.null(names(data))) dataset_names = paste0("Dataset_", 1:length(data)) else dataset_names = names(data)
      Ns = sapply(data, ncol)
      data_list = data
      
    }
    
  }
  
  if (datatype == 'incidence_raw') {
    
    if (is.null(names(data))) dataset_names = paste0("Dataset_", 1:length(data)) else dataset_names = names(data)
    Ns = sapply(data, length)
    data_list = data
    
  }
  
  
  if (is.null(conf)) conf = 0.95
  tmp = qnorm(1 - (1 - conf)/2)
  
  trunc = ifelse(is.null(level), T, F)
  
  if (base == 'coverage' ) {
    
    if ( is.null(level) ) {
      
      level <- lapply(1:length(data_list), function(i) {
        
        level = seq(0.5, 1, 0.025)
        
        if(datatype == "abundance") {
          
          n = sum(data_list[[i]])
          data_gamma = rowSums(data_list[[i]])
          data_gamma = data_gamma[data_gamma>0]
          data_alpha = as.matrix(data_list[[i]]) %>% as.vector
          
          ref_gamma = iNEXT.3D:::Coverage(data_gamma, 'abundance', n)
          ref_alpha = iNEXT.3D:::Coverage(data_alpha, 'abundance', n)
          ref_alpha_max = iNEXT.3D:::Coverage(data_alpha, 'abundance', n*2)
          ref_gamma_max = iNEXT.3D:::Coverage(data_gamma, 'abundance', n*2)
          
        }else if(datatype == "incidence_raw"){
          
          n = unique(sapply(data_list[[i]], ncol))
          gamma = Reduce('+', data_list[[i]])
          gamma[gamma>1] = 1
          data_gamma_raw = gamma
          data_gamma_freq = c(n, rowSums(gamma))
          
          data_alpha_freq = sapply(data_list[[i]], rowSums) %>% c(n, .)
          
          ref_gamma = iNEXT.3D:::Coverage(data_gamma_freq, 'incidence_freq', n)
          ref_alpha = iNEXT.3D:::Coverage(data_alpha_freq, 'incidence_freq', n)
          ref_alpha_max = iNEXT.3D:::Coverage(data_alpha_freq, 'incidence_freq', n*2)
          ref_gamma_max = iNEXT.3D:::Coverage(data_gamma_freq, 'incidence_freq', n*2)
          
        }
        
        c(level, ref_gamma, ref_alpha, ref_alpha_max, ref_gamma_max) %>% sort %>% unique
        
      })
      
    } else {
      
      if (class(level) == "numeric" | class(level) == "integer" | class(level) == "double") {
        level <- list(level = level)
      }
      
      if (length(level) != length(data_list)) level <- lapply(1:length(data_list), function(x) level[[1]])
      
    }
    
    } else if ( base == 'size' ) {
    
      if ( is.null(level) ) {
        
        if (datatype == "abundance") {
          
          endpoint <- sapply(data_list, function(x) 2*sum(x))
          
        } else if (datatype == "incidence_raw") {
          
          endpoint <- sapply(data_list, function(x) 2*ncol(x[[1]]))
          
        }
        
        level <- lapply(1:length(data_list), function(i) {
          
          if(datatype == "abundance") {
            
            ni <- sum(data_list[[i]])
            
          }else if(datatype == "incidence_raw"){
            
            ni <- ncol(data_list[[i]][[1]])
          }
          
          mi <- floor(c(seq(1, ni-1, length.out = 20), ni, seq(ni+1, endpoint[i], length.out = 20)))
        })
        
      } else {
        
        if (class(level) == "numeric" | class(level) == "integer" | class(level) == "double") {
          level <- list(level = level)
        }
        
        if (length(level) != length(data_list)) level <- lapply(1:length(data_list), function(x) level[[1]])
        
      }
  }
  
  if (diversity == 'FD' & FDtype == 'tau_value' & is.null(FDtau) == T) {
    
    if (datatype == 'abundance') {
      
      tt <- lapply(data_list, rowSums)
      tt = lapply(tt, function(i) data.frame('value' = i) %>% rownames_to_column(var = "Species"))
      pdata = tt[[1]]
      
      if (length(tt) > 1) {
        for(i in 2:length(tt)){
          pdata = full_join(pdata, tt[[i]], by = "Species")
        }
      }
      
      pdata[is.na(pdata)] = 0
      pdata = pdata %>% column_to_rownames("Species")
      pdata = rowSums(pdata)
      
      order_sp <- match(names(pdata),rownames(FDdistM))
      FDdistM <- FDdistM[order_sp,order_sp]
      pdata <- matrix(pdata/sum(pdata), ncol = 1)
      
    } else if (datatype == 'incidence_raw') {
      
      tt <- lapply(data_list, function(x) {tt = Reduce('+', x); tt[tt > 1] = 1; rowSums(tt) })
      tt = lapply(tt, function(i) data.frame('value' = i) %>% rownames_to_column(var = "Species"))
      pdata = tt[[1]]
      
      if (length(tt) > 1) {
        for(i in 2:length(tt)){
          pdata = full_join(pdata, tt[[i]], by = "Species")
        }
      }
      
      pdata[is.na(pdata)] = 0
      pdata = pdata %>% column_to_rownames("Species")
      pdata = rowSums(pdata)
      
      order_sp <- match(names(pdata),rownames(FDdistM))
      FDdistM <- FDdistM[order_sp,order_sp]
      pdata <- matrix(pdata/sum(pdata), ncol = 1)
      
    }
    
    FDtau <- sum ( (pdata %*% t(pdata) ) * FDdistM) # dmean
  }
  
  if (diversity == 'PD') {
    
    if (datatype == "abundance") 
      
      if (length(data_list) > 1) {
        
        pool.data = data_list[[1]] %>% data.frame %>% rownames_to_column()
        for (i in 2:length(data_list)) 
          pool.data = full_join(pool.data, data_list[[i]] %>% data.frame %>% rownames_to_column(), 'rowname')
        pool.data[is.na(pool.data)] = 0
        pool.data = pool.data %>% column_to_rownames() %>% rowSums
        
      } else pool.data = do.call(cbind, data_list) %>% rowSums
    
    if (datatype == 'incidence_raw') pool.data = do.call(cbind,lapply(data_list, function(x) do.call(cbind,x)) ) %>% rowSums
    
    pool.name = names(pool.data[pool.data>0])
    tip = PDtree$tip.label[-match(pool.name, PDtree$tip.label)]
    mytree = drop.tip(PDtree, tip)
    H_max = get.rooted.tree.height(mytree)
    
    if(is.null(PDreftime)) { reft = H_max
    } else if (PDreftime <= 0) { stop("Reference time must be greater than 0. Use NULL to set it to pooled tree height.", call. = FALSE)
    } else { reft = PDreftime }
    
  }
    
  for_each_dataset = function(data, dataset_name, N, level) {
    
    #data
    if (datatype == 'abundance') {
      
      n = sum(data)
      data_gamma = rowSums(data)
      data_gamma = data_gamma[data_gamma>0]
      data_alpha = as.matrix(data) %>% as.vector
      
      ref_gamma = iNEXT.3D:::Coverage(data_gamma, 'abundance', n)
      ref_alpha = iNEXT.3D:::Coverage(data_alpha, 'abundance', n)
      ref_alpha_max = iNEXT.3D:::Coverage(data_alpha, 'abundance', n*2)
      ref_gamma_max = iNEXT.3D:::Coverage(data_gamma, 'abundance', n*2)
      
      # level = c(level, ref_gamma, ref_alpha, ref_alpha_max, ref_gamma_max) %>% sort %>% unique
      # level = level[level<1]
      
      m_gamma = sapply(level, function(i) coverage_to_size(data_gamma, i, datatype='abundance'))
      m_alpha = sapply(level, function(i) coverage_to_size(data_alpha, i, datatype='abundance'))
      
    }
    
    if (datatype == 'incidence_raw') {
      
      sampling_units = sapply(data, ncol)
      if (length(unique(sampling_units)) > 1) stop("unsupported data structure: the sampling units of all datasets must be the same.")
      if (length(unique(sampling_units)) == 1) n = unique(sampling_units)
      
      gamma = Reduce('+', data)
      gamma[gamma>1] = 1
      data_gamma_raw = gamma
      data_gamma_freq = c(n, rowSums(gamma))
      
      data_alpha_freq = sapply(data, rowSums) %>% c(n, .)
      
      # data_gamma_freq = data_gamma_freq[data_gamma_freq>0]
      # data_alpha_freq = data_alpha_freq[data_alpha_freq>0]
      
      data_2D = apply(sapply(data, rowSums), 2, function(x) c(n, x)) %>% as.data.frame
      
      ref_gamma = iNEXT.3D:::Coverage(data_gamma_freq, 'incidence_freq', n)
      ref_alpha = iNEXT.3D:::Coverage(data_alpha_freq, 'incidence_freq', n)
      ref_alpha_max = iNEXT.3D:::Coverage(data_alpha_freq, 'incidence_freq', n*2)
      ref_gamma_max = iNEXT.3D:::Coverage(data_gamma_freq, 'incidence_freq', n*2)
      
      # level = c(level, ref_gamma, ref_alpha, ref_alpha_max, ref_gamma_max) %>% sort %>% unique
      # level = level[level < 1]
      
      m_gamma = sapply(level, function(i) coverage_to_size(data_gamma_freq, i, datatype='incidence_freq'))
      m_alpha = sapply(level, function(i) coverage_to_size(data_alpha_freq, i, datatype='incidence_freq'))
      
      # if (sum(m_gamma < 1 | m_alpha < 1) > 0) {
      #   index = m_gamma < 1 | m_alpha < 1
      #   level = level[!index]
      #   m_gamma = m_gamma[!index]
      #   m_alpha = m_alpha[!index]
      # }
    }
    
    
    
    if (diversity == 'TD') {
      
      if (datatype == 'abundance') {
        
        gamma = lapply(1:length(level), function(i){
          estimate3D(as.numeric(data_gamma), diversity = 'TD', q = q, datatype = "abundance", base = "coverage", level = level[i], nboot = 0)
        }) %>% do.call(rbind,.)
        
        alpha = lapply(1:length(level), function(i){
          estimate3D(as.numeric(data_alpha), diversity = 'TD', q = q, datatype = "abundance", base = "coverage", level = level[i], nboot = 0)
        }) %>% do.call(rbind,.)
        
      }
      
      if (datatype == 'incidence_raw') {
        
        gamma = lapply(1:length(level), function(i){
          estimate3D(as.numeric(data_gamma_freq), diversity = 'TD', q = q, datatype = "incidence_freq", base = "coverage", level = level[i], nboot = 0)
        }) %>% do.call(rbind,.)
        
        alpha = lapply(1:length(level), function(i){
          estimate3D(data_alpha_freq, diversity = 'TD', q = q, datatype = "incidence_freq", base = "coverage", level = level[i], nboot = 0)
        }) %>% do.call(rbind,.)
        
        
        
      }
      
      gamma = (cbind(SC = rep(level, each=length(q)), gamma[,-c(1,3,7,8,9)]) %>% 
                 mutate(Method = ifelse(SC>=ref_gamma, ifelse(SC == ref_gamma, 'Observed', 'Extrapolation'), 'Rarefaction'))
      )[,c(4,2,5,1,3)] %>% set_colnames(c('Estimate', 'Order.q', 'Method', 'SC', 'Size'))
      
      
      if (max_alpha_coverage == T) under_max_alpha = !((gamma$Order.q == 0) & (gamma$SC > ref_alpha_max)) else under_max_alpha = gamma$SC > 0
      gamma = gamma[under_max_alpha,]
      
      
      
      alpha = (cbind(SC = rep(level, each = length(q)), alpha[,-c(1,3,7,8,9)]) %>% 
                 mutate(Method = ifelse(SC >= ref_alpha, ifelse(SC == ref_alpha, 'Observed', 'Extrapolation'), 'Rarefaction'))
      )[,c(4,2,5,1,3)] %>% set_colnames(c('Estimate', 'Order.q', 'Method', 'SC', 'Size'))
      
      alpha$Estimate = alpha$Estimate / N
      
      
      alpha = alpha[under_max_alpha,]
      
      
      
      beta = alpha
      beta$Estimate = gamma$Estimate/alpha$Estimate
      # beta[beta == "Observed"] = "Observed_SC(n, alpha)"
      # beta = beta %>% rbind(., cbind(gamma %>% filter(Method == "Observed") %>% select(Estimate) / alpha %>% filter(Method == "Observed") %>% select(Estimate), 
      #                                Order.q = q, Method = "Observed", SC = NA, Size = beta[beta$Method == "Observed_SC(n, alpha)", 'Size']))
      
      C = beta %>% mutate(Estimate = ifelse(Order.q==1, log(Estimate)/log(N), (Estimate^(1-Order.q) - 1)/(N^(1-Order.q)-1)))
      U = beta %>% mutate(Estimate = ifelse(Order.q==1, log(Estimate)/log(N), (Estimate^(Order.q-1) - 1)/(N^(Order.q-1)-1)))
      V = beta %>% mutate(Estimate = (Estimate-1)/(N-1))
      S = beta %>% mutate(Estimate = (1/Estimate-1)/(1/N-1))
      
      if(nboot>1){
        
        # cl = makeCluster(cluster_numbers)
        # clusterExport(cl, c("bootstrap_population_multiple_assemblage","data","data_gamma", 'data_gamma_freq',"level","N",'under_max_alpha',
        #                     'datatype', 'data_2D'))
        # clusterEvalQ(cl, library(tidyverse, magrittr))
        
        # plan(sequential)
        # plan(multiprocess)
        
        # se = parSapply(cl, 1:nboot, function(i){
        
        # start = Sys.time()
        se = future_lapply(1:nboot, function(i){
          
          if (datatype == 'abundance') {
            
            bootstrap_population = bootstrap_population_multiple_assemblage(data, data_gamma, 'abundance')
            bootstrap_sample = sapply(1:ncol(data), function(k) rmultinom(n = 1, size = sum(data[,k]), prob = bootstrap_population[,k]))
            
            bootstrap_data_gamma = rowSums(bootstrap_sample)
            bootstrap_data_gamma = bootstrap_data_gamma[bootstrap_data_gamma > 0]
            bootstrap_data_alpha = as.matrix(bootstrap_sample) %>% as.vector
            bootstrap_data_alpha = bootstrap_data_alpha[bootstrap_data_alpha > 0]
            
            gamma = lapply(1:length(level), function(i){
              estimate3D(as.numeric(bootstrap_data_gamma), diversity = 'TD', q = q, datatype = "abundance", base = "coverage", level = level[i], nboot = 0)
            }) %>% do.call(rbind,.)
            
            alpha = lapply(1:length(level), function(i){
              estimate3D(as.numeric(bootstrap_data_alpha), diversity = 'TD', q = q, datatype = "abundance", base = "coverage", level = level[i], nboot = 0)
            }) %>% do.call(rbind,.)
            
          }
          
          if (datatype == 'incidence_raw') {
            
            bootstrap_population = bootstrap_population_multiple_assemblage(data_2D, data_gamma_freq, 'incidence')
            
            raw = lapply(1:ncol(bootstrap_population), function(j){
              
              lapply(1:nrow(bootstrap_population), function(i) rbinom(n = n, size = 1, prob = bootstrap_population[i,j])) %>% do.call(rbind,.)
              
            })
            
            gamma = Reduce('+', raw)
            gamma[gamma > 1] = 1
            bootstrap_data_gamma_freq = c(n, rowSums(gamma))
            
            bootstrap_data_alpha_freq = sapply(raw, rowSums) %>% c(n, .)
            
            bootstrap_data_gamma_freq = bootstrap_data_gamma_freq[bootstrap_data_gamma_freq > 0]
            bootstrap_data_alpha_freq = bootstrap_data_alpha_freq[bootstrap_data_alpha_freq > 0]
            
            gamma = lapply(1:length(level), function(i){
              estimate3D(bootstrap_data_gamma_freq, diversity = 'TD', q = q, datatype = "incidence_freq", base = "coverage", level = level[i], nboot = 0)
            }) %>% do.call(rbind,.)
            
            alpha = lapply(1:length(level), function(i){
              estimate3D(bootstrap_data_alpha_freq, diversity = 'TD', q = q, datatype = "incidence_freq", base = "coverage", level = level[i], nboot = 0)
            }) %>% do.call(rbind,.)
            
          }
          
          gamma = gamma$qTD[under_max_alpha]
          
          alpha = alpha$qTD[under_max_alpha]
          alpha = alpha / N
          
          beta = gamma/alpha
          
          Order.q = rep(q, length(level))[under_max_alpha]
          
          beta = data.frame(Estimate=beta, Order.q)
          
          C = (beta %>% mutate(Estimate = ifelse(Order.q == 1,log(Estimate)/log(N),(Estimate^(1 - Order.q) - 1)/(N^(1 - Order.q) - 1))))$Estimate
          U = (beta %>% mutate(Estimate = ifelse(Order.q == 1,log(Estimate)/log(N),(Estimate^(Order.q - 1) - 1)/(N^(Order.q - 1) - 1))))$Estimate
          V = (beta %>% mutate(Estimate = (Estimate - 1)/(N - 1)))$Estimate
          S = (beta %>% mutate(Estimate = (1/Estimate - 1)/(1/N - 1)))$Estimate
          
          beta = beta$Estimate
          
          cbind(gamma, alpha, beta, C, U, V, S) %>% as.matrix
          
          # }, simplify = "array") %>% apply(., 1:2, sd) %>% data.frame
        }, future.seed = TRUE) %>% abind(along = 3) %>% apply(1:2, sd)
        # end = Sys.time()
        # end - start
        
        # stopCluster(cl)
        # plan(sequential)
        
      } else {
        
        se = matrix(NA, ncol = 7, nrow = nrow(beta))
        colnames(se) = c("gamma", "alpha", "beta", "C", "U", 'V', 'S')
        se = as.data.frame(se)
        
      }
      
      if (nboot>0 & datatype == 'abundance') {
        gamma.se = estimate3D(data_gamma, diversity = 'TD', q = q, datatype = datatype, base = "coverage", level = level, nboot = nboot) %>% arrange(., SC, Order.q) %>% select(s.e.) %>% unlist
        
        alpha.se = estimate3D(data_alpha, diversity = 'TD', q = q, datatype = datatype, base = "coverage", level = level, nboot = nboot) %>% arrange(., SC, Order.q) %>% select(s.e.) %>% unlist
        alpha.se = alpha.se / N
        
      } else if (nboot>0 & datatype == 'incidence_raw') {
        
        gamma.se = estimate3D(data_gamma_freq, diversity = 'TD', q = q, datatype = 'incidence_freq', base = "coverage", level = level, nboot = nboot) %>% arrange(., SC, Order.q) %>% select(s.e.) %>% unlist
        
        alpha.se = estimate3D(data_alpha_freq, diversity = 'TD', q = q, datatype = 'incidence_freq', base = "coverage", level = level, nboot = nboot) %>% arrange(., SC, Order.q) %>% select(s.e.) %>% unlist
        alpha.se = alpha.se / N
        
      } else {
        gamma.se = alpha.se = NA
      }
      
      se[1:( length(level) * length(q) ), 'gamma'] = gamma.se
      
      se[1:( length(level) * length(q) ), 'alpha'] = alpha.se
      
    }
    
    if (diversity == 'PD') {
      
      if (datatype == 'abundance') {
        
        aL = iNEXT.3D:::phyBranchAL_Abu(phylo = PDtree, data = data_gamma, rootExtend = T, refT = reft)
        aL$treeNabu$branch.length = aL$BLbyT[,1]
        aL_table_gamma = aL$treeNabu %>% select(branch.abun, branch.length, tgroup)
        
        gamma = iNEXT.3D:::PhD.m.est(ai = aL_table_gamma$branch.abun, Lis = as.matrix(aL_table_gamma$branch.length), m = m_gamma, q = q, nt = n, cal = "PD") %>% t %>% as.data.frame %>% 
          set_colnames(q) %>% gather(Order.q, Estimate) %>% 
          mutate(SC = rep(level, length(q)), Size = rep(m_gamma, length(q)))
        
        
        aL_table_alpha = c()
        
        for (i in 1:N){
          
          x = data[data[,i]>0,i]
          names(x) = rownames(data)[data[,i]>0]
          
          aL = iNEXT.3D:::phyBranchAL_Abu(phylo = PDtree, data = x, rootExtend = T, refT = reft)
          aL$treeNabu$branch.length = aL$BLbyT[,1]
          aL_table = aL$treeNabu %>% select(branch.abun, branch.length, tgroup)
          
          aL_table_alpha = rbind(aL_table_alpha, aL_table)
          
        }
        
        
        qPDm = iNEXT.3D:::PhD.m.est(ai = aL_table_alpha$branch.abun, Lis = as.matrix(aL_table_alpha$branch.length), m = m_alpha, q = q, nt = n, cal = "PD")
        qPDm = qPDm/N
        alpha = qPDm %>% t %>% as.data.frame %>% 
          set_colnames(q) %>% gather(Order.q, Estimate) %>% 
          mutate(SC = rep(level, length(q)), Size = rep(m_alpha, length(q)))
        
      }
      
      if (datatype == 'incidence_raw') {
        
        aL = iNEXT.3D:::phyBranchAL_Inc(phylo = PDtree, data = as.matrix(data_gamma_raw), datatype = "incidence_raw", refT = reft, rootExtend = T)
        aL$treeNabu$branch.length = aL$BLbyT[,1]
        aL_table_gamma = aL$treeNabu %>% select(branch.abun, branch.length, tgroup)
        gamma = iNEXT.3D:::PhD.m.est(ai = aL_table_gamma$branch.abun, Lis=as.matrix(aL_table_gamma$branch.length), m = m_gamma, q = q, nt = n, cal = "PD") %>% t %>% as.data.frame %>% 
          set_colnames(q) %>% gather(Order.q, Estimate) %>% 
          mutate(SC = rep(level, length(q)), Size = rep(m_gamma, length(q)))
        
        aL_table_alpha = c()
        
        for (i in 1:N){
          
          x = data[[i]]
          
          aL = iNEXT.3D:::phyBranchAL_Inc(phylo = PDtree, data = x, datatype = "incidence_raw", rootExtend = T, refT = reft)
          aL$treeNabu$branch.length = aL$BLbyT[,1]
          aL_table = aL$treeNabu %>% select(branch.abun, branch.length, tgroup)
          
          aL_table_alpha = rbind(aL_table_alpha, aL_table)
          
        }
        
        alpha = (iNEXT.3D:::PhD.m.est(ai = aL_table_alpha$branch.abun, Lis = as.matrix(aL_table_alpha$branch.length), m = m_alpha, q = q, nt = n, cal = "PD")/N) %>% t %>% as.data.frame %>% 
          set_colnames(q) %>% gather(Order.q, Estimate) %>% 
          mutate(SC = rep(level, length(q)), Size = rep(m_alpha, length(q)))
        
        
      }
      
      gamma = (gamma %>% 
                 mutate(Method = ifelse(SC >= ref_gamma, ifelse(SC == ref_gamma, 'Observed', 'Extrapolation'), 'Rarefaction')))[,c(2,1,5,3,4)] %>% 
        set_colnames(c('Estimate', 'Order.q', 'Method', 'SC', 'Size'))
      
      if (max_alpha_coverage == T) under_max_alpha = !((gamma$Order.q == 0) & (gamma$SC > ref_alpha_max)) else under_max_alpha = gamma$SC>0
      gamma = gamma[under_max_alpha,]
      gamma$Order.q = as.numeric(gamma$Order.q)
      
      
      alpha = (alpha %>% 
                 mutate(Method = ifelse(SC >= ref_alpha, ifelse(SC == ref_alpha, 'Observed', 'Extrapolation'), 'Rarefaction')))[,c(2,1,5,3,4)] %>% 
        set_colnames(c('Estimate', 'Order.q', 'Method', 'SC', 'Size'))
      
      alpha = alpha[under_max_alpha,]
      alpha$Order.q = as.numeric(alpha$Order.q)
      
      if (PDtype == 'meanPD') {
        gamma$Estimate = gamma$Estimate/reft
        alpha$Estimate = alpha$Estimate/reft
      }
      
      beta = alpha
      beta$Estimate = gamma$Estimate/alpha$Estimate
      # beta[beta == "Observed"] = "Observed_SC(n, alpha)"
      # beta = beta %>% rbind(., cbind(gamma %>% filter(Method == "Observed") %>% select(Estimate) / alpha %>% filter(Method == "Observed") %>% select(Estimate), 
      #                                Order.q = q, Method = "Observed", SC = NA, Size = beta[beta$Method == "Observed_SC(n, alpha)", 'Size']))
      
      C = beta %>% mutate(Estimate = ifelse(Order.q == 1, log(Estimate)/log(N), (Estimate^(1 - Order.q) - 1)/(N^(1 - Order.q) - 1)))
      U = beta %>% mutate(Estimate = ifelse(Order.q == 1, log(Estimate)/log(N), (Estimate^(Order.q - 1) - 1)/(N^(Order.q - 1) - 1)))
      V = beta %>% mutate(Estimate = (Estimate - 1)/(N - 1))
      S = beta %>% mutate(Estimate = (1/Estimate - 1)/(1/N - 1))
      
      if(nboot>1){
        
        # cl = makeCluster(cluster_numbers)
        # clusterExport(cl, c("bootstrap_population_multiple_assemblage","data","data_gamma", 'data_gamma_freq',"level","N",'under_max_alpha',
        #                     'datatype', 'data_2D'))
        # clusterEvalQ(cl, library(tidyverse, magrittr))
        
        # plan(sequential)
        # plan(multiprocess)
        
        # se = parSapply(cl, 1:nboot, function(i){
        
        # start = Sys.time()
        se = future_lapply(1:nboot, function(i){
          
          if (datatype == 'abundance') {
            
            tree_bt = PDtree
            
            bootstrap_population = bootstrap_population_multiple_assemblage(data, data_gamma, 'abundance')
            p_bt = bootstrap_population
            unseen_p = p_bt[-(1:nrow(data)),] %>% matrix(ncol = ncol(data))
            
            if ( nrow(p_bt) > nrow(data) & sum(unseen_p) > 0 ){
              
              unseen = unseen_p[which(rowSums(unseen_p) > 0),]
              unseen = matrix(unseen, ncol = ncol(unseen_p), byrow = T)
              p_bt = rbind(p_bt[(1:nrow(data)),], unseen)
              unseen_name = sapply(1:nrow(unseen), function(i) paste0('unseen_', i))
              rownames(p_bt) = c(rownames(data), unseen_name)
              
              bootstrap_sample = sapply(1:ncol(data), function(k) rmultinom(n = 1, size = sum(data[,k]), prob = p_bt[,k]))
              x_bt = bootstrap_sample
              
              rownames(x_bt) = rownames(p_bt)
              
              if ( sum(x_bt[-(1:nrow(data)),])>0 ){
                
                g0_hat = apply(data, 2, function(x){
                  
                  n = sum(x)
                  f1 = sum(x == 1)
                  f2 = sum(x == 2)
                  
                  aL = iNEXT.3D:::phyBranchAL_Abu(phylo = PDtree, data = x, rootExtend = T, refT = reft)
                  
                  aL$treeNabu$branch.length = aL$BLbyT[,1]
                  aL = aL$treeNabu %>% select(branch.abun,branch.length)
                  g1 = aL$branch.length[aL$branch.abun == 1] %>% sum
                  g2 = aL$branch.length[aL$branch.abun == 2] %>% sum
                  g0_hat = ifelse( g2 > ((g1*f2)/(2*f1)) , ((n-1)/n)*(g1^2/(2*g2)) , ((n-1)/n)*(g1*(f1-1)/(2*(f2+1))) )
                  if(is.na(g0_hat)) {g0_hat <- 0 }
                  g0_hat
                  
                })
                
                te = (x_bt[1:nrow(data),]*(data == 0))>0
                used_length = sapply(1:ncol(data), function(i) { 
                  
                  if (sum(te[,i]) == 0) return(0) else {
                    
                    iNEXT.3D:::phyBranchAL_Abu(phylo = PDtree, data = x_bt[1:nrow(data),i], rootExtend = T, refT = reft)$treeNabu %>%
                      subset(label %in% names(which(te[,i] == TRUE))) %>% select(branch.length) %>% sum
                    
                  }
                  
                })
                
                g0_hat = g0_hat - used_length
                g0_hat[g0_hat < 0] = 0
                
                unseen_sample = x_bt[-(1:nrow(data)),]
                if (is.vector(unseen_sample)) unseen_sample = matrix(unseen_sample, ncol = ncol(x_bt), byrow = T)
                
                L0_hat = sapply(1:length(g0_hat), function(i) if(sum(unseen_sample[,i] > 0) > 0) (g0_hat[i] / nrow(unseen)) else 0 )
                
                L0_hat = rowSums((matrix(L0_hat, nrow(unseen_sample), ncol(unseen_sample), byrow = T) * unseen_sample)) / rowSums(unseen_sample)
                L0_hat[which(rowSums(unseen_sample) == 0)] = 0
                
                for (i in 1:length(L0_hat)){
                  
                  tip = list(edge = matrix(c(2,1),1,2),
                             tip.label = unseen_name[i],
                             edge.length = L0_hat[i],
                             Nnode = 1)
                  class(tip) = "phylo"
                  
                  tree_bt = tree_bt + tip
                  
                }
                
              } else {
                
                x_bt = x_bt[1:nrow(data),]
                p_bt = p_bt[1:nrow(data),]
                
              }
              
            } else {
              
              p_bt = p_bt[1:nrow(data),]
              x_bt = sapply(1:ncol(data), function(k) rmultinom(n = 1, size = sum(data[,k]), prob = p_bt[,k]))
              rownames(x_bt) = rownames(data)
              
            }
            
            bootstrap_data_gamma = rowSums(x_bt)
            bootstrap_data_gamma = bootstrap_data_gamma[bootstrap_data_gamma>0]
            bootstrap_data_alpha = as.matrix(x_bt) %>% as.vector
            bootstrap_data_alpha = bootstrap_data_alpha[bootstrap_data_alpha>0]
            
            m_gamma = sapply(level, function(i) coverage_to_size(bootstrap_data_gamma, i, datatype='abundance'))
            m_alpha = sapply(level, function(i) coverage_to_size(bootstrap_data_alpha, i, datatype='abundance'))
            
            aL = iNEXT.3D:::phyBranchAL_Abu(phylo = tree_bt, data = bootstrap_data_gamma, rootExtend = T, refT = reft)
            aL$treeNabu$branch.length = aL$BLbyT[,1]
            aL_table_gamma = aL$treeNabu %>% select(branch.abun, branch.length, tgroup)
            
            gamma = as.vector(iNEXT.3D:::PhD.m.est(ai = aL_table_gamma$branch.abun, Lis = as.matrix(aL_table_gamma$branch.length), m = m_gamma, q = q, nt = n, cal = "PD") %>% t)
            
            
            aL_table_alpha = c()
            
            for (i in 1:N){
              
              # x = x_bt[x_bt[,i]>0,i]
              # names(x) = rownames(p_bt)[x_bt[,i]>0]
              
              x = x_bt[,i]
              names(x) = rownames(p_bt)
              x = x[x_bt[,i]>0]
              
              aL = iNEXT.3D:::phyBranchAL_Abu(phylo = tree_bt, data = x, rootExtend = T, refT = reft)
              aL$treeNabu$branch.length = aL$BLbyT[,1]
              aL_table = aL$treeNabu %>% select(branch.abun, branch.length, tgroup)
              
              aL_table_alpha = rbind(aL_table_alpha, aL_table)
              
            }
            
            alpha = as.vector((iNEXT.3D:::PhD.m.est(ai = aL_table_alpha$branch.abun, Lis = as.matrix(aL_table_alpha$branch.length), m = m_alpha, q = q, nt = n, cal = "PD")/N) %>% t)
            
            # beta_obs = (iNEXT.3D:::PD.Tprofile(ai = aL_table_gamma$branch.abun, Lis = as.matrix(aL_table_gamma$branch.length), q = q, nt = n, cal = "PD") / 
            #               (iNEXT.3D:::PD.Tprofile(ai = aL_table_alpha$branch.abun, Lis = as.matrix(aL_table_alpha$branch.length), q = q, nt = n, cal = "PD") / N)) %>% unlist()
          }
          
          if (datatype == 'incidence_raw') {
            
            tree_bt = PDtree
            
            bootstrap_population = bootstrap_population_multiple_assemblage(data_2D, data_gamma_freq, 'incidence')
            p_bt = bootstrap_population
            unseen_p = p_bt[-(1:nrow(data[[1]])),] %>% matrix(ncol=N)
            
            if ( nrow(p_bt) > nrow(data[[1]]) & sum(unseen_p)>0 ){
              
              unseen = unseen_p[which(rowSums(unseen_p)>0),]
              unseen = matrix(unseen, ncol = ncol(unseen_p), byrow = T)
              p_bt = rbind(p_bt[(1:nrow(data[[1]])),], unseen)
              unseen_name = sapply(1:nrow(unseen), function(i) paste0('unseen_', i))
              rownames(p_bt) = c(rownames(data[[1]]), unseen_name)
              
              raw = lapply(1:ncol(p_bt), function(j){
                
                lapply(1:nrow(p_bt), function(i) rbinom(n=n, size=1, prob=p_bt[i,j])) %>% do.call(rbind,.)
                
              })
              
              for (i in 1:length(raw)) rownames(raw[[i]]) = rownames(p_bt)
              
              if ( lapply(1:length(raw), function(i) raw[[i]][-(1:nrow(data[[1]])),]) %>% do.call(sum,.)>0 ){
                
                R0_hat = sapply(data, function(x){
                  
                  nT = ncol(x)
                  Q1 = sum(rowSums(x)==1)
                  Q2 = sum(rowSums(x)==2)
                  
                  aL = iNEXT.3D:::phyBranchAL_Inc(phylo = PDtree, data = x, datatype = "incidence_raw", rootExtend = T, refT = reft)
                  
                  aL$treeNabu$branch.length = aL$BLbyT[,1]
                  aL = aL$treeNabu %>% select(branch.abun,branch.length)
                  R1 = aL$branch.length[aL$branch.abun == 1] %>% sum
                  R2 = aL$branch.length[aL$branch.abun == 2] %>% sum
                  R0_hat = ifelse( R2>((R1*Q2)/(2*Q1)) , ((nT-1)/nT)*(R1^2/(2*R2)) , ((nT-1)/nT)*(R1*(Q1-1)/(2*(Q2+1))) )
                  if(is.na(R0_hat)) { R0_hat <- 0 }
                  R0_hat
                  
                })
                
                te = (sapply(raw, rowSums)[1:nrow(data[[1]]),]*(sapply(data, rowSums) == 0)) > 0
                used_length = sapply(1:N, function(i) {
                  
                  if (sum(te[,i]) == 0) return(0) else {
                    
                    iNEXT.3D:::phyBranchAL_Inc(phylo = PDtree, data = raw[[i]][1:nrow(data[[1]]),], datatype = "incidence_raw", rootExtend = T, refT = reft)$treeNabu %>%
                      subset(label %in% names(which(te[,i] == TRUE))) %>% select(branch.length) %>% sum
                    
                  }
                  
                })
                
                R0_hat = R0_hat - used_length
                R0_hat[R0_hat < 0] = 0
                
                unseen_sample = sapply(raw, rowSums)[-(1:nrow(data[[1]])),]
                if (is.vector(unseen_sample)) unseen_sample = matrix(unseen_sample, ncol = N, byrow = T)
                
                L0_hat = sapply(1:length(R0_hat), function(i) if(sum(unseen_sample[,i] > 0) > 0) (R0_hat[i] / nrow(unseen)) else 0 )
                
                L0_hat = rowSums((matrix(L0_hat, nrow(unseen_sample), ncol(unseen_sample), byrow = T) * unseen_sample)) / rowSums(unseen_sample)
                L0_hat[which(rowSums(unseen_sample) == 0)] = 0
                
                for (i in 1:length(L0_hat)){
                  
                  tip = list(edge = matrix(c(2,1), 1, 2),
                             tip.label = unseen_name[i],
                             edge.length = L0_hat[i],
                             Nnode = 1)
                  class(tip) = "phylo"
                  
                  tree_bt = tree_bt + tip
                  
                }
                
              } else raw = lapply(raw, function(i) i[1:nrow(data[[1]]),])
              
            } else {
              
              p_bt = p_bt[1:nrow(data[[1]]),]
              raw = lapply(1:ncol(p_bt), function(j){
                
                lapply(1:nrow(p_bt), function(i) rbinom(n = n, size = 1, prob = p_bt[i,j])) %>% do.call(rbind,.)
                
              })
              
              for (i in 1:length(raw)) rownames(raw[[i]]) = rownames(p_bt)
              
            }
            
            gamma = Reduce('+', raw)
            gamma[gamma > 1] = 1
            bootstrap_data_gamma_raw = gamma
            bootstrap_data_gamma_freq = c(n, rowSums(gamma))
            
            bootstrap_data_alpha_freq = sapply(raw, rowSums) %>% c(n, .)
            
            bootstrap_data_gamma_freq = bootstrap_data_gamma_freq[bootstrap_data_gamma_freq > 0]
            bootstrap_data_alpha_freq = bootstrap_data_alpha_freq[bootstrap_data_alpha_freq > 0]
            
            m_gamma = sapply(level, function(i) coverage_to_size(bootstrap_data_gamma_freq, i, datatype = 'incidence'))
            m_alpha = sapply(level, function(i) coverage_to_size(bootstrap_data_alpha_freq, i, datatype = 'incidence'))
            
            aL = iNEXT.3D:::phyBranchAL_Inc(phylo = tree_bt, data = bootstrap_data_gamma_raw, datatype = "incidence_raw", rootExtend = T, refT = reft)
            aL$treeNabu$branch.length = aL$BLbyT[,1]
            aL_table_gamma = aL$treeNabu %>% select(branch.abun, branch.length, tgroup)
            
            gamma = as.vector(iNEXT.3D:::PhD.m.est(ai = aL_table_gamma$branch.abun, Lis = as.matrix(aL_table_gamma$branch.length), m = m_gamma, q = q, nt = n, cal="PD") %>% t)
            
            
            aL_table_alpha = c()
            
            for (i in 1:N){
              
              x = raw[[i]]
              
              aL = iNEXT.3D:::phyBranchAL_Inc(phylo = tree_bt, data = x, datatype = "incidence_raw", rootExtend = T, refT = reft)
              aL$treeNabu$branch.length = aL$BLbyT[,1]
              aL_table = aL$treeNabu %>% select(branch.abun, branch.length, tgroup)
              
              aL_table_alpha = rbind(aL_table_alpha, aL_table)
              
            }
            
            alpha = as.vector((iNEXT.3D:::PhD.m.est(ai = aL_table_alpha$branch.abun, Lis = as.matrix(aL_table_alpha$branch.length), m = m_alpha, q=q, nt = n, cal = "PD")/N) %>% t)
            
            # beta_obs = (iNEXT.3D:::PD.Tprofile(ai = aL_table_gamma$branch.abun, Lis = as.matrix(aL_table_gamma$branch.length), q = q, nt = n, cal = "PD") / 
            #               (iNEXT.3D:::PD.Tprofile(ai = aL_table_alpha$branch.abun, Lis = as.matrix(aL_table_alpha$branch.length), q = q, nt = n, cal = "PD") / N)) %>% unlist()
          }
          
          gamma = gamma[under_max_alpha]
          
          alpha = alpha[under_max_alpha]
          
          if (PDtype == 'meanPD') {
            gamma = gamma/reft
            alpha = alpha/reft
          } 
          
          beta = gamma/alpha
          
          # gamma = c(gamma, rep(0, length(q)))
          # alpha = c(alpha, rep(0, length(q)))
          # beta = c(beta, beta_obs)
          # 
          # Order.q = rep(q, each = length(level) + 1)[under_max_alpha]
          Order.q = rep(q, each = length(level))[under_max_alpha]
          
          beta = data.frame(Estimate = beta, Order.q)
          
          C = (beta %>% mutate(Estimate = ifelse(Order.q == 1, log(Estimate)/log(N), (Estimate^(1 - Order.q) - 1)/(N^(1 - Order.q) - 1))))$Estimate
          U = (beta %>% mutate(Estimate = ifelse(Order.q == 1, log(Estimate)/log(N), (Estimate^(Order.q - 1) - 1)/(N^(Order.q - 1) - 1))))$Estimate
          V = (beta %>% mutate(Estimate = (Estimate - 1)/(N - 1)))$Estimate
          S = (beta %>% mutate(Estimate = (1/Estimate - 1)/(1/N - 1)))$Estimate
          
          beta = beta$Estimate
          
          cbind(gamma, alpha, beta, C, U, V, S) %>% as.matrix
          
          # }, simplify = "array") %>% apply(., 1:2, sd) %>% data.frame
        }, future.seed = TRUE) %>% abind(along = 3) %>% apply(1:2, sd)
        # end = Sys.time()
        # end - start
        
        # stopCluster(cl) 
        # plan(sequential)
        
      } else {
        
        se = matrix(NA, ncol = 7, nrow = nrow(beta))
        colnames(se) = c("gamma", "alpha", "beta", "C", "U", 'V', 'S')
        se = as.data.frame(se)
        
      }
      
    }
    
    if (diversity == 'FD') {
      
      FDdistM = as.matrix(FDdistM)
      
      FD_by_tau = function(data, distM, tau, level, datatype, by, m_gamma, m_alpha) {
        
        if (datatype == 'abundance') {
          
          zik = data
          zik = zik[rowSums(data)>0,]
          
          dij = distM
          dij = dij[rowSums(data)>0, rowSums(data)>0]
          
          
          if (tau == 0) {
            dij[dij > 0] <- 1
            aik = (1 - dij/1) %*% as.matrix(zik)
          } else {
            dij[which(dij>tau, arr.ind = T)] = tau
            aik = (1 - dij/tau) %*% as.matrix(zik)
          }
          positive_id = rowSums(aik)>0
          
          
          gamma_x = rowSums(zik)[positive_id]
          gamma_a = rowSums(aik)[positive_id]
          gamma_v = gamma_x/gamma_a
          # gamma_a = ifelse(gamma_a < 1, 1, round(gamma_a))
          
          ai_vi_gamma = list(ai = data.frame(gamma_a), vi = data.frame(gamma_v))
          
          gamma = FD.m.est_0(ai_vi_gamma, m_gamma, q, n) %>% as.vector
          # gamma = iNEXT.3D:::FD.m.est(ai_vi_gamma, m_gamma, q, n) %>% as.vector
          
          
          alpha_x = as.vector(as.matrix(zik))
          alpha_a = as.vector(aik)
          # alpha_a = ifelse(alpha_a < 1, 1, round(alpha_a))
          
          alpha_v = alpha_x/alpha_a
          alpha_v = rep(gamma_v,N)
          
          alpha_v = alpha_v[alpha_a>0]
          alpha_a = alpha_a[alpha_a>0]
          
          ai_vi_alpha = list(ai = data.frame(alpha_a), vi = data.frame(alpha_v))
          
          if (by == 'size') alpha = (FD.m.est_0(ai_vi_alpha, m_gamma, q, n)/N) %>% as.vector
          if (by == 'coverage') alpha = (FD.m.est_0(ai_vi_alpha, m_alpha, q, n)/N) %>% as.vector
          # if (by == 'size') alpha = (iNEXT.3D:::FD.m.est(ai_vi_alpha, m_gamma, q, n)/N) %>% as.vector
          # if (by == 'coverage') alpha = (iNEXT.3D:::FD.m.est(ai_vi_alpha, m_alpha, q, n)/N) %>% as.vector
          
        }
        
        if (datatype == 'incidence_raw') {
          
          data_gamma_freq = data$data_gamma_freq
          data_2D = data$data_2D
          
          gamma_Y = data_gamma_freq[-1]
          
          dij = distM
          dij = dij[gamma_Y > 0, gamma_Y > 0]
          gamma_Y = gamma_Y[gamma_Y > 0]
          
          
          if (tau == 0) {
            dij[dij > 0] <- 1
            gamma_a = (1 - dij/1) %*% as.matrix(gamma_Y)
          } else {
            dij[which(dij>tau, arr.ind = T)] = tau
            gamma_a = (1 - dij/tau) %*% as.matrix(gamma_Y)
          }
          
          gamma_a[gamma_a > n] = n
          gamma_v = gamma_Y/gamma_a
          
          # gamma_a = ifelse(gamma_a < 1, 1, round(gamma_a))
          # gamma_a[gamma_a > n] = n
          
          ai_vi_gamma = list(ai = data.frame(gamma_a), vi = data.frame(gamma_v))
          
          gamma = FD.m.est_0(ai_vi_gamma, m_gamma, q, n) %>% as.vector
          # gamma = iNEXT.3D:::FD.m.est(ai_vi_gamma, m_gamma, q, n) %>% as.vector
          
          
          alpha_Y = data_2D[-1,]
          
          dij = distM
          dij = dij[rowSums(data_2D[-1,]) > 0, rowSums(data_2D[-1,])>0]
          alpha_Y = alpha_Y[rowSums(data_2D[-1,])>0,]
          
          
          if (tau == 0) {
            dij[dij > 0] <- 1
            alpha_a = (1 - dij/1) %*% as.matrix(alpha_Y)
          } else {
            dij[which(dij>tau, arr.ind = T)] = tau
            alpha_a = (1 - dij/tau) %*% as.matrix(alpha_Y)
          }
          
          # alpha_a = ifelse(alpha_a < 1, 1, round(alpha_a))
          alpha_a[alpha_a > n] = n
          alpha_a = as.vector(alpha_a)
          
          alpha_v = rep(gamma_v, N)
          alpha_v = alpha_v[alpha_a > 0]
          alpha_a = alpha_a[alpha_a > 0]
          
          ai_vi_alpha = list(ai = data.frame(alpha_a), vi = data.frame(alpha_v))
          
          if (by == 'size') alpha = (FD.m.est_0(ai_vi_alpha, m_gamma, q, n)/N) %>% as.vector
          if (by == 'coverage') alpha = (FD.m.est_0(ai_vi_alpha, m_alpha, q, n)/N) %>% as.vector
          # if (by == 'size') alpha = (iNEXT.3D:::FD.m.est(ai_vi_alpha, m_gamma, q, n)/N) %>% as.vector
          # if (by == 'coverage') alpha = (iNEXT.3D:::FD.m.est(ai_vi_alpha, m_alpha, q, n)/N) %>% as.vector
          
        }
        
        return(data.frame(gamma,alpha))
        
      }
      
      if (FDtype == 'tau_value'){
        
        if (datatype == 'abundance') {
          
          FDdistM = FDdistM[rownames(FDdistM) %in% rownames(data), colnames(FDdistM) %in% rownames(data)]
          order_sp <- match(rownames(data),rownames(FDdistM))
          FDdistM <- FDdistM[order_sp,order_sp]
          
          output = FD_by_tau(data, FDdistM, FDtau, level, datatype = 'abundance', by = 'coverage', m_gamma = m_gamma, m_alpha = m_alpha)
          gamma = output$gamma
          alpha = output$alpha
          
          gamma = data.frame(SC = level, gamma) %>% 
            mutate(Method = ifelse(SC>=ref_gamma, ifelse(SC == ref_gamma, 'Observed', 'Extrapolation'), 'Rarefaction'),
                   Order.q = rep(q, each=length(SC)/length(q)), Size = rep(m_gamma, length(q)))
          
          alpha = data.frame(SC = level, alpha) %>% 
            mutate(Method = ifelse(SC>=ref_alpha, ifelse(SC == ref_alpha, 'Observed', 'Extrapolation'), 'Rarefaction'),
                   Order.q = rep(q, each=length(SC)/length(q)), Size = rep(m_alpha, length(q)))
          
          beta = alpha
          beta$alpha = gamma$gamma/alpha$alpha
          
          ## Observed Beta ##
          # output_obs = FD_by_tau(data, FDdistM, FDtau, level, datatype = 'abundance', by = 'size', m_gamma = sum(data), m_alpha = sum(data))
          # gamma_obs = output_obs$gamma
          # alpha_obs = output_obs$alpha
          # obs_beta = gamma_obs/alpha_obs
          
        }
        
        if (datatype == 'incidence_raw') {
          
          FDdistM = FDdistM[rownames(FDdistM) %in% names(data_gamma_freq)[-1], colnames(FDdistM) %in% names(data_gamma_freq)[-1]]
          order_sp <- match(names(data_gamma_freq)[-1],rownames(FDdistM))
          FDdistM <- FDdistM[order_sp,order_sp]
          
          output = FD_by_tau(list(data_gamma_freq = data_gamma_freq, data_2D = data_2D), FDdistM, FDtau, level, datatype='incidence_raw', by = 'coverage', m_gamma=m_gamma, m_alpha=m_alpha)
          gamma = output$gamma
          alpha = output$alpha
          
          gamma = data.frame(SC = level, gamma) %>% 
            mutate(Method = ifelse(SC >= ref_gamma, ifelse(SC == ref_gamma, 'Observed', 'Extrapolation'), 'Rarefaction'),
                   Order.q = rep(q, each = length(SC)/length(q)), Size = rep(m_gamma, length(q)))
          
          alpha = data.frame(SC = level, alpha) %>% 
            mutate(Method = ifelse(SC>=ref_alpha, ifelse(SC == ref_alpha, 'Observed', 'Extrapolation'), 'Rarefaction'),
                   Order.q = rep(q, each = length(SC)/length(q)), Size = rep(m_alpha, length(q)))
          
          beta = alpha
          beta$alpha = gamma$gamma/alpha$alpha
          
          ## Observed Beta ##
          # output_obs = FD_by_tau(list(data_gamma_freq = data_gamma_freq, data_2D = data_2D), FDdistM, FDtau, level, datatype = 'incidence_raw', by = 'size', m_gamma = data_gamma_freq[1], m_alpha = data_gamma_freq[1])
          # gamma_obs = output_obs$gamma
          # alpha_obs = output_obs$alpha
          # obs_beta = gamma_obs/alpha_obs
          
        }
        
      }
      
      if (FDtype == 'AUC'){
        
        # cut = seq(0.00000001, 1, length.out = FDcut_number)
        cut = seq(0, 1, length.out = FDcut_number)
        width = diff(cut)
        
        if (datatype == 'abundance') {
          
          FDdistM = FDdistM[rownames(FDdistM) %in% rownames(data), colnames(FDdistM) %in% rownames(data)]
          order_sp <- match(rownames(data),rownames(FDdistM))
          FDdistM <- FDdistM[order_sp,order_sp]
          
          gamma_alpha_over_tau = lapply(cut, function(tau) {
            
            FD_by_tau(data, FDdistM, tau, level, datatype = 'abundance', by = 'coverage', m_gamma = m_gamma, m_alpha = m_alpha)
            
          })
          
          gamma_over_tau = sapply(gamma_alpha_over_tau, function(x) x$gamma)
          
          left_limit  = apply(gamma_over_tau, 1, function(x) x[-FDcut_number]*width)
          right_limit = apply(gamma_over_tau, 1, function(x) x[-1]*width)
          
          gamma = colSums((left_limit + right_limit)/2)
          
          alpha_over_tau = sapply(gamma_alpha_over_tau, function(x) x$alpha)
          
          left_limit  = apply(alpha_over_tau, 1, function(x) x[-FDcut_number]*width)
          right_limit = apply(alpha_over_tau, 1, function(x) x[-1]*width)
          
          alpha = colSums((left_limit + right_limit)/2)
          
          beta_over_tau = gamma_over_tau/alpha_over_tau
          
          left_limit  = apply(beta_over_tau, 1, function(x) x[-FDcut_number]*width)
          right_limit = apply(beta_over_tau, 1, function(x) x[-1]*width)
          
          beta = colSums((left_limit + right_limit)/2)
          
          gamma = data.frame(SC = level, gamma) %>% 
            mutate(Method = ifelse(SC >= ref_gamma, ifelse(SC == ref_gamma, 'Observed', 'Extrapolation'), 'Rarefaction'),
                   Order.q = rep(q, each = length(SC)/length(q)), Size = rep(m_gamma, length(q)))
          
          alpha = data.frame(SC = level, alpha) %>% 
            mutate(Method = ifelse(SC >= ref_alpha, ifelse(SC == ref_alpha, 'Observed', 'Extrapolation'), 'Rarefaction'),
                   Order.q = rep(q, each = length(SC)/length(q)), Size = rep(m_alpha, length(q)))
          
          beta = data.frame(SC = level, beta) %>% 
            mutate(Method = ifelse(SC >= ref_alpha, ifelse(SC == ref_alpha, 'Observed', 'Extrapolation'), 'Rarefaction'),
                   Order.q = rep(q, each = length(SC)/length(q)), Size = rep(m_alpha, length(q)))
          
          ## Observed Beta ##
          # obs_gamma_alpha_over_tau = lapply(cut, function(tau) {
          #   
          #   FD_by_tau(data, FDdistM, tau, level, datatype = 'abundance', by = 'size', m_gamma = sum(data), m_alpha = sum(data))
          #   
          # })
          # 
          # obs_beta_over_tau = sapply(obs_gamma_alpha_over_tau, function(x) x$gamma) / sapply(obs_gamma_alpha_over_tau, function(x) x$alpha)
          # 
          # if (length(q) == 1) obs_beta_over_tau = matrix(obs_beta_over_tau, nrow = 1)
          # 
          # obs_beta = colSums( (apply(obs_beta_over_tau, 1, function(x) x[-FDcut_number]*width) + apply(obs_beta_over_tau, 1, function(x) x[-1]*width) ) / 2)
          
        }
        
        if (datatype == 'incidence_raw') {
          
          FDdistM = FDdistM[rownames(FDdistM) %in% names(data_gamma_freq)[-1], colnames(FDdistM) %in% names(data_gamma_freq)[-1]]
          order_sp <- match(names(data_gamma_freq)[-1],rownames(FDdistM))
          FDdistM <- FDdistM[order_sp,order_sp]
          
          gamma_alpha_over_tau = lapply(cut, function(tau) {
            
            FD_by_tau(list(data_gamma_freq = data_gamma_freq, data_2D = data_2D), FDdistM, tau, level, datatype = 'incidence_raw', by = 'coverage', m_gamma = m_gamma, m_alpha = m_alpha)
            
          })
          
          gamma_over_tau = sapply(gamma_alpha_over_tau, function(x) x$gamma)
          
          left_limit  = apply(gamma_over_tau, 1, function(x) x[-FDcut_number]*width)
          right_limit = apply(gamma_over_tau, 1, function(x) x[-1]*width)
          
          gamma = colSums((left_limit + right_limit)/2)
          
          alpha_over_tau = sapply(gamma_alpha_over_tau, function(x) x$alpha)
          
          left_limit  = apply(alpha_over_tau, 1, function(x) x[-FDcut_number]*width)
          right_limit = apply(alpha_over_tau, 1, function(x) x[-1]*width)
          
          alpha = colSums((left_limit + right_limit)/2)
          
          beta_over_tau = gamma_over_tau/alpha_over_tau
          
          left_limit  = apply(beta_over_tau, 1, function(x) x[-FDcut_number]*width)
          right_limit = apply(beta_over_tau, 1, function(x) x[-1]*width)
          
          beta = colSums((left_limit + right_limit)/2)
          
          gamma = data.frame(SC = level, gamma) %>% 
            mutate(Method = ifelse(SC >= ref_gamma, ifelse(SC == ref_gamma, 'Observed', 'Extrapolation'), 'Rarefaction'),
                   Order.q = rep(q, each = length(SC)/length(q)), Size = rep(m_gamma, length(q)))
          
          alpha = data.frame(SC = level, alpha) %>% 
            mutate(Method = ifelse(SC >= ref_alpha, ifelse(SC == ref_alpha, 'Observed', 'Extrapolation'), 'Rarefaction'),
                   Order.q = rep(q, each = length(SC)/length(q)), Size = rep(m_alpha, length(q)))
          
          beta = data.frame(SC = level, beta) %>% 
            mutate(Method = ifelse(SC >= ref_alpha, ifelse(SC == ref_alpha, 'Observed', 'Extrapolation'), 'Rarefaction'),
                   Order.q = rep(q, each = length(SC)/length(q)), Size = rep(m_alpha, length(q)))
          
          ## Observed Beta ##
          # obs_gamma_alpha_over_tau = lapply(cut, function(tau) {
          #   
          #   FD_by_tau(list(data_gamma_freq = data_gamma_freq, data_2D = data_2D), FDdistM, tau, level, datatype = 'incidence_raw', by = 'size', m_gamma = data_gamma_freq[1], m_alpha = data_alpha_freq[1])
          #   
          # })
          # 
          # obs_beta_over_tau = sapply(obs_gamma_alpha_over_tau, function(x) x$gamma) / sapply(obs_gamma_alpha_over_tau, function(x) x$alpha)
          # 
          # if (length(q) == 1) obs_beta_over_tau = matrix(obs_beta_over_tau, nrow = 1)
          # 
          # obs_beta = colSums( (apply(obs_beta_over_tau, 1, function(x) x[-FDcut_number]*width) + apply(obs_beta_over_tau, 1, function(x) x[-1]*width) ) / 2)
          
        }
        
      }
      
      gamma = gamma[,c(2,4,3,1,5)] %>% set_colnames(c('Estimate', 'Order.q', 'Method', 'SC', 'Size'))
      
      if (max_alpha_coverage == T) under_max_alpha = !((gamma$Order.q == 0) & (gamma$SC > ref_alpha_max)) else under_max_alpha = gamma$SC > 0
      gamma = gamma[under_max_alpha,]
      
      
      alpha = alpha[,c(2,4,3,1,5)] %>% set_colnames(c('Estimate', 'Order.q', 'Method', 'SC', 'Size'))
      
      alpha = alpha[under_max_alpha,]
      
      beta = beta[,c(2,4,3,1,5)] %>% set_colnames(c('Estimate', 'Order.q', 'Method', 'SC', 'Size'))
      
      beta = beta[under_max_alpha,]
      
      # beta[beta == "Observed"] = "Observed_SC(n, alpha)"
      # beta = beta %>% 
      #   rbind(., data.frame(Estimate = obs_beta, Order.q = q, Method = "Observed", SC = NA, Size = beta[beta$Method == "Observed_SC(n, alpha)", 'Size']))
      
      C = beta %>% mutate(Estimate = ifelse(Order.q == 1, log(Estimate)/log(N), (Estimate^(1 - Order.q) - 1)/(N^(1 - Order.q) - 1)))
      U = beta %>% mutate(Estimate = ifelse(Order.q == 1, log(Estimate)/log(N), (Estimate^(Order.q - 1) - 1)/(N^(Order.q - 1) - 1)))
      V = beta %>% mutate(Estimate = (Estimate - 1)/(N - 1))
      S = beta %>% mutate(Estimate = (1/Estimate - 1)/(1/N - 1))
      
      if(nboot > 1){
        
        # cl = makeCluster(cluster_numbers)
        # clusterExport(cl, c("bootstrap_population_multiple_assemblage","data","data_gamma", 'data_gamma_freq',"level","N",'under_max_alpha',
        #                     'datatype', 'data_2D'))
        # clusterEvalQ(cl, library(tidyverse, magrittr))
        
        # plan(sequential)
        # plan(multiprocess, workers=7)
        
        # se = parSapply(cl, 1:nboot, function(i){
        
        # start = Sys.time()
        se = future_lapply(1:nboot, function(i){
          
          if (datatype == 'abundance') {
            
            p_bt = bootstrap_population_multiple_assemblage(data, data_gamma, 'abundance')
            f0_hat = nrow(p_bt) - nrow(data)
            
            distance_matrix_bt = Bootstrap_distance_matrix(rowSums(data), FDdistM, f0_hat, 'abundance')
            
            data_bt = sapply(1:ncol(data), function(k) rmultinom(n = 1, size = sum(data[,k]), prob = p_bt[,k]))
            
            data_gamma = rowSums(data_bt)
            data_gamma = data_gamma[data_gamma>0]
            data_alpha = as.matrix(data_bt) %>% as.vector
            
            if (FDtype == 'tau_value'){
              
              m_gamma = sapply(level, function(i) coverage_to_size(data_gamma, i, datatype='abundance'))
              m_alpha = sapply(level, function(i) coverage_to_size(data_alpha, i, datatype='abundance'))
              
              output = FD_by_tau(data_bt, distance_matrix_bt, FDtau, level, datatype='abundance', by = 'coverage', m_gamma = m_gamma, m_alpha = m_alpha)
              gamma = output$gamma
              alpha = output$alpha
              
              beta=gamma/alpha
              
              ## Observed Beta ##
              # output_obs = FD_by_tau(data_bt, distance_matrix_bt, FDtau, level, datatype = 'abundance', by = 'size', m_gamma = sum(data_bt), m_alpha = sum(data_bt))
              # gamma_obs = output_obs$gamma
              # alpha_obs = output_obs$alpha
              # beta_obs = gamma_obs/alpha_obs
              
            }
            
            if (FDtype == 'AUC'){
              
              m_gamma = sapply(level, function(i) coverage_to_size(data_gamma, i, datatype='abundance'))
              m_alpha = sapply(level, function(i) coverage_to_size(data_alpha, i, datatype='abundance'))
              
              gamma_alpha_over_tau = lapply(cut, function(tau) {
                
                FD_by_tau(data_bt, distance_matrix_bt, tau, level, datatype = 'abundance', by = 'coverage', m_gamma = m_gamma, m_alpha = m_alpha)
                
              })
              
              gamma_over_tau = sapply(gamma_alpha_over_tau, function(x) x$gamma)
              
              left_limit  = apply(gamma_over_tau, 1, function(x) x[-FDcut_number]*width)
              right_limit = apply(gamma_over_tau, 1, function(x) x[-1]*width)
              
              gamma = colSums((left_limit + right_limit)/2)
              
              alpha_over_tau = sapply(gamma_alpha_over_tau, function(x) x$alpha)
              
              left_limit  = apply(alpha_over_tau, 1, function(x) x[-FDcut_number]*width)
              right_limit = apply(alpha_over_tau, 1, function(x) x[-1]*width)
              
              alpha = colSums((left_limit + right_limit)/2)
              
              beta_over_tau = gamma_over_tau/alpha_over_tau
              
              left_limit  = apply(beta_over_tau, 1, function(x) x[-FDcut_number]*width)
              right_limit = apply(beta_over_tau, 1, function(x) x[-1]*width)
              
              beta = colSums((left_limit + right_limit)/2)
              
              ## Observed Beta ##
              # gamma_alpha_over_tau = lapply(cut, function(tau) {
              #   
              #   FD_by_tau(data_bt, distance_matrix_bt, tau, level, datatype = 'abundance', by = 'size', m_gamma = sum(data_bt), m_alpha = sum(data_bt))
              #   
              # })
              # 
              # gamma_over_tau = sapply(gamma_alpha_over_tau, function(x) x$gamma)
              # 
              # alpha_over_tau = sapply(gamma_alpha_over_tau, function(x) x$alpha)
              # 
              # beta_over_tau = gamma_over_tau/alpha_over_tau
              # 
              # if (length(q) == 1) beta_over_tau = matrix(beta_over_tau, nrow = 1)
              # 
              # left_limit  = apply(beta_over_tau, 1, function(x) x[-FDcut_number]*width)
              # right_limit = apply(beta_over_tau, 1, function(x) x[-1]*width)
              # 
              # beta_obs = colSums((left_limit + right_limit)/2)
              
            }
            
          }
          
          if (datatype == 'incidence_raw') {
            
            p_bt = bootstrap_population_multiple_assemblage(data_2D, data_gamma_freq, 'incidence')
            f0_hat = nrow(p_bt) - nrow(data_2D[-1,])
            
            distance_matrix_bt = Bootstrap_distance_matrix(c(n,rowSums(data_gamma_raw)), FDdistM, f0_hat, 'incidence_freq')
            
            # p_bt = p_bt[rowSums(p_bt)>0,]
            
            raw = lapply(1:ncol(p_bt), function(j){
              
              lapply(1:nrow(p_bt), function(i) rbinom(n = n, size = 1, prob = p_bt[i,j])) %>% do.call(rbind,.)
              
            })
            
            gamma = Reduce('+', raw)
            gamma[gamma>1] = 1
            data_gamma_raw_bt = gamma
            data_gamma_freq_bt = c(n, rowSums(gamma))
            
            data_alpha_freq_bt = sapply(raw, rowSums) %>% c(n, .)
            
            # data_gamma_freq_bt = data_gamma_freq_bt[data_gamma_freq_bt > 0]
            # data_alpha_freq_bt = data_alpha_freq_bt[data_alpha_freq_bt > 0]
            
            data_2D_bt = apply(sapply(raw, rowSums), 2, function(x) c(n, x)) %>% as.data.frame
            
            if (FDtype == 'tau_value'){
              
              m_gamma = sapply(level, function(i) coverage_to_size(data_gamma_freq_bt, i, datatype='incidence_freq'))
              m_alpha = sapply(level, function(i) coverage_to_size(data_alpha_freq_bt, i, datatype='incidence_raw'))
              
              output = FD_by_tau(list(data_gamma_freq = data_gamma_freq_bt, data_2D = data_2D_bt), distance_matrix_bt, FDtau, level, datatype = 'incidence_raw', by = 'coverage', m_gamma = m_gamma, m_alpha = m_alpha)
              gamma = output$gamma
              alpha = output$alpha
              
              beta = gamma/alpha
              
              ## Observed Beta ##
              # output_obs = FD_by_tau(list(data_gamma_freq = data_gamma_freq_bt, data_2D = data_2D_bt), distance_matrix_bt, FDtau, level, datatype = 'incidence_raw', by = 'size', m_gamma = data_gamma_freq_bt[1], m_alpha = data_gamma_freq_bt[1])
              # gamma_obs = output_obs$gamma
              # alpha_obs = output_obs$alpha
              # beta_obs = gamma_obs/alpha_obs
              
            }
            
            if (FDtype == 'AUC'){
              
              m_gamma = sapply(level, function(i) coverage_to_size(data_gamma_freq_bt, i, datatype='incidence_freq'))
              m_alpha = sapply(level, function(i) coverage_to_size(data_alpha_freq_bt, i, datatype='incidence_raw'))
              
              gamma_alpha_over_tau = lapply(cut, function(tau) {
                
                FD_by_tau(list(data_gamma_freq = data_gamma_freq_bt, data_2D = data_2D_bt), distance_matrix_bt, tau, level, datatype = 'incidence_raw', by = 'coverage', m_gamma = m_gamma, m_alpha = m_alpha)
                
              })
              
              gamma_over_tau = sapply(gamma_alpha_over_tau, function(x) x$gamma)
              
              left_limit  = apply(gamma_over_tau, 1, function(x) x[-FDcut_number]*width)
              right_limit = apply(gamma_over_tau, 1, function(x) x[-1]*width)
              
              gamma = colSums((left_limit + right_limit)/2)
              
              alpha_over_tau = sapply(gamma_alpha_over_tau, function(x) x$alpha)
              
              left_limit  = apply(alpha_over_tau, 1, function(x) x[-FDcut_number]*width)
              right_limit = apply(alpha_over_tau, 1, function(x) x[-1]*width)
              
              alpha = colSums((left_limit + right_limit)/2)
              
              beta_over_tau = gamma_over_tau/alpha_over_tau
              
              left_limit  = apply(beta_over_tau, 1, function(x) x[-FDcut_number]*width)
              right_limit = apply(beta_over_tau, 1, function(x) x[-1]*width)
              
              beta = colSums((left_limit + right_limit)/2)
              
              ## Observed Beta ##
              # gamma_alpha_over_tau = lapply(cut, function(tau) {
              #   
              #   FD_by_tau(list(data_gamma_freq = data_gamma_freq_bt, data_2D = data_2D_bt), distance_matrix_bt, tau, level, datatype = 'incidence_raw', by = 'size', m_gamma = data_gamma_freq_bt[1], m_alpha = data_gamma_freq_bt[1])
              #   
              # })
              # 
              # gamma_over_tau = sapply(gamma_alpha_over_tau, function(x) x$gamma)
              # 
              # alpha_over_tau = sapply(gamma_alpha_over_tau, function(x) x$alpha)
              # 
              # beta_over_tau = gamma_over_tau/alpha_over_tau
              # 
              # if (length(q) == 1) beta_over_tau = matrix(beta_over_tau, nrow = 1)
              # 
              # left_limit  = apply(beta_over_tau, 1, function(x) x[-FDcut_number]*width)
              # right_limit = apply(beta_over_tau, 1, function(x) x[-1]*width)
              # 
              # beta_obs = colSums((left_limit + right_limit)/2)
              
            }
            
          }
          
          gamma = gamma[under_max_alpha]
          alpha = alpha[under_max_alpha]
          beta = beta[under_max_alpha]
          
          beta = gamma/alpha
          
          # gamma = c(gamma, rep(0, length(q)))
          # alpha = c(alpha, rep(0, length(q)))
          # beta = c(beta, beta_obs)
          # 
          # Order.q = rep(q, each=length(level) + 1)[under_max_alpha]
          Order.q = rep(q, each=length(level))[under_max_alpha]
          
          beta = data.frame(Estimate = beta, Order.q)
          
          C = (beta %>% mutate(Estimate = ifelse(Order.q == 1, log(Estimate)/log(N), (Estimate^(1 - Order.q) - 1)/(N^(1 - Order.q) - 1))))$Estimate
          U = (beta %>% mutate(Estimate = ifelse(Order.q == 1, log(Estimate)/log(N), (Estimate^(Order.q - 1) - 1)/(N^(Order.q - 1) - 1))))$Estimate
          V = (beta %>% mutate(Estimate = (Estimate - 1)/(N - 1)))$Estimate
          S = (beta %>% mutate(Estimate = (1/Estimate - 1)/(1/N - 1)))$Estimate
          
          beta = beta$Estimate
          
          cbind(gamma, alpha, beta, C, U, V, S) %>% as.matrix
          
          # }, simplify = "array") %>% apply(., 1:2, sd) %>% data.frame
        }, future.seed = TRUE) %>% abind(along = 3) %>% apply(1:2, sd)
        # end = Sys.time()
        # end - start
        
        # stopCluster(cl)
        # plan(sequential)
        
      } else {
        
        se = matrix(NA, ncol = 7, nrow = nrow(beta))
        colnames(se) = c("gamma", "alpha", "beta", "C", "U", 'V', 'S')
        se = as.data.frame(se)
        
      }
      
    }
    
    
    se = as.data.frame(se)
    
    if (diversity == "TD") index = "TD"
    if (diversity == "PD" & PDtype == "PD") index = "PD"
    if (diversity == "PD" & PDtype == "meanPD") index = "meanPD"
    if (diversity == "FD" & FDtype == "tau_value") index = "FD_tau"
    if (diversity == "FD" & FDtype == "AUC") index = "FD_AUC"
    
    gamma = gamma %>% mutate(s.e. = se$gamma,
                             LCL = Estimate - tmp * se$gamma,
                             UCL = Estimate + tmp * se$gamma,
                             Dataset = dataset_name,
                             Diversity = index) %>% 
      arrange(Order.q, SC) %>% .[,c(9, 2, 4, 5, 1, 3, 6, 7, 8, 10)] %>% rename("Gamma" = "Estimate")
    
    alpha = alpha %>% mutate(s.e. = se$alpha,
                             LCL = Estimate - tmp * se$alpha,
                             UCL = Estimate + tmp * se$alpha,
                             Dataset = dataset_name,
                             Diversity = index) %>% 
      arrange(Order.q, SC) %>% .[,c(9, 2, 4, 5, 1, 3, 6, 7, 8, 10)] %>% rename("Alpha" = "Estimate")
    
    beta = beta %>% mutate(  s.e. = se$beta,
                             LCL = Estimate - tmp * se$beta,
                             UCL = Estimate + tmp * se$beta,
                             Dataset = dataset_name,
                             Diversity = index) %>% 
      arrange(Order.q, SC) %>% .[,c(9, 2, 4, 5, 1, 3, 6, 7, 8, 10)] %>% rename("Beta" = "Estimate")
    
    C = C %>% mutate(        s.e. = se$C,
                             LCL = Estimate - tmp * se$C,
                             UCL = Estimate + tmp * se$C,
                             Dataset = dataset_name,
                             Diversity = index) %>% 
      arrange(Order.q, SC) %>% .[,c(9, 2, 4, 5, 1, 3, 6, 7, 8, 10)] %>% rename("Dissimilarity" = "Estimate")
    
    
    U = U %>% mutate(        s.e. = se$U,
                             LCL = Estimate - tmp * se$U,
                             UCL = Estimate + tmp * se$U,
                             Dataset = dataset_name,
                             Diversity = index) %>% 
      arrange(Order.q, SC) %>% .[,c(9, 2, 4, 5, 1, 3, 6, 7, 8, 10)] %>% rename("Dissimilarity" = "Estimate")
    
    V = V %>% mutate(        s.e. = se$V,
                             LCL = Estimate - tmp * se$V,
                             UCL = Estimate + tmp * se$V,
                             Dataset = dataset_name,
                             Diversity = index) %>% 
      arrange(Order.q, SC) %>% .[,c(9, 2, 4, 5, 1, 3, 6, 7, 8, 10)] %>% rename("Dissimilarity" = "Estimate")
    
    S = S %>% mutate(        s.e. = se$S,
                             LCL = Estimate - tmp * se$S,
                             UCL = Estimate + tmp * se$S,
                             Dataset = dataset_name,
                             Diversity = index) %>% 
      arrange(Order.q, SC) %>% .[,c(9, 2, 4, 5, 1, 3, 6, 7, 8, 10)] %>% rename("Dissimilarity" = "Estimate")
    
    
    if (trunc) {
      
      gamma = gamma %>% filter(!(Order.q==0 & round(Size)>2*n))
      
      alpha = alpha %>% filter(!(Order.q==0 & round(Size)>2*n))
      
      beta  = beta  %>% filter(!(Order.q==0 & round(Size)>2*n))
      
      C    =  C    %>% filter(!(Order.q==0 & round(Size)>2*n))
      
      U    =  U    %>% filter(!(Order.q==0 & round(Size)>2*n))
      
      V    =  V    %>% filter(!(Order.q==0 & round(Size)>2*n))
      
      S    =  S    %>% filter(!(Order.q==0 & round(Size)>2*n))
      
    }
    
    
    if (diversity == "FD" & FDtype == "tau_value") {
      
      gamma = gamma %>% mutate(Tau = FDtau)
      
      alpha = alpha %>% mutate(Tau = FDtau)
      
      beta  = beta  %>% mutate(Tau = FDtau)
      
      C    =  C    %>% mutate(Tau = FDtau)
      
      U    =  U    %>% mutate(Tau = FDtau)
      
      V    =  V    %>% mutate(Tau = FDtau)
      
      S    =  S    %>% mutate(Tau = FDtau)
      
    }
    
    
    if (datatype == "abundance") {
      
      alpha$Method[alpha$SC == ref_gamma_max] = 
        beta$Method[beta$SC == ref_gamma_max] = 
        C$Method[C$SC == ref_gamma_max] = 
        U$Method[U$SC == ref_gamma_max] = 
        V$Method[V$SC == ref_gamma_max] = 
        S$Method[S$SC == ref_gamma_max] = "Extrap_SC(2n, gamma)"
      
      gamma$Method[gamma$SC == ref_alpha_max] = 
        alpha$Method[alpha$SC == ref_alpha_max] = 
        beta$Method[beta$SC == ref_alpha_max] = 
        C$Method[C$SC == ref_alpha_max] = 
        U$Method[U$SC == ref_alpha_max] = 
        V$Method[V$SC == ref_alpha_max] = 
        S$Method[S$SC == ref_alpha_max] = "Extrap_SC(2n, alpha)"
      
      gamma$Method[gamma$SC == ref_gamma_max] = "Extrap_SC(2n, gamma)"
      
      
      alpha$Method[alpha$SC == ref_gamma] = 
        beta$Method[beta$SC == ref_gamma] = 
        C$Method[C$SC == ref_gamma] = 
        U$Method[U$SC == ref_gamma] = 
        V$Method[V$SC == ref_gamma] = 
        S$Method[S$SC == ref_gamma] = "Observed_SC(n, gamma)"
      
      gamma$Method[gamma$SC == ref_alpha] = 
        alpha$Method[alpha$SC == ref_alpha] = 
        beta$Method[beta$SC == ref_alpha] = 
        C$Method[C$SC == ref_alpha] = 
        U$Method[U$SC == ref_alpha] = 
        V$Method[V$SC == ref_alpha] = 
        S$Method[S$SC == ref_alpha] = "Observed_SC(n, alpha)"
      
      gamma$Method[gamma$SC == ref_gamma] = "Observed_SC(n, gamma)"
      
      
    }
    
    if (datatype == "incidence_raw") {
      
      alpha$Method[alpha$SC == ref_gamma_max] = 
        beta$Method[beta$SC == ref_gamma_max] = 
        C$Method[C$SC == ref_gamma_max] = 
        U$Method[U$SC == ref_gamma_max] = 
        V$Method[V$SC == ref_gamma_max] = 
        S$Method[S$SC == ref_gamma_max] = "Extrap_SC(2T, gamma)"
      
      gamma$Method[gamma$SC == ref_alpha_max] = 
        alpha$Method[alpha$SC == ref_alpha_max] = 
        beta$Method[beta$SC == ref_alpha_max] = 
        C$Method[C$SC == ref_alpha_max] = 
        U$Method[U$SC == ref_alpha_max] = 
        V$Method[V$SC == ref_alpha_max] = 
        S$Method[S$SC == ref_alpha_max] = "Extrap_SC(2T, alpha)"
      
      gamma$Method[gamma$SC == ref_gamma_max] = "Extrap_SC(2T, gamma)"
      
      
      alpha$Method[alpha$SC == ref_gamma] = 
        beta$Method[beta$SC == ref_gamma] = 
        C$Method[C$SC == ref_gamma] = 
        U$Method[U$SC == ref_gamma] = 
        V$Method[V$SC == ref_gamma] = 
        S$Method[S$SC == ref_gamma] = "Observed_SC(T, gamma)"
      
      gamma$Method[gamma$SC == ref_alpha] = 
        alpha$Method[alpha$SC == ref_alpha] = 
        beta$Method[beta$SC == ref_alpha] = 
        C$Method[C$SC == ref_alpha] = 
        U$Method[U$SC == ref_alpha] = 
        V$Method[V$SC == ref_alpha] = 
        S$Method[S$SC == ref_alpha] = "Observed_SC(T, alpha)"
      
      gamma$Method[gamma$SC == ref_gamma] = "Observed_SC(T, gamma)"
      
    }
    
    
    
    list(gamma = gamma, alpha = alpha, beta = beta, `1-C` = C, `1-U` = U, `1-V` = V, `1-S` = S)
    
  }
  
  for_each_dataset.size = function(data, dataset_name, N, level) {
    
    #data
    if (datatype == 'abundance') {
      
      n = sum(data)
      data_gamma = rowSums(data)
      data_gamma = data_gamma[data_gamma>0]
      data_alpha = as.matrix(data) %>% as.vector
      
      ref_gamma = n
      ref_alpha = n
      
    }
    
    if (datatype == 'incidence_raw') {
      
      sampling_units = sapply(data, ncol)
      if (length(unique(sampling_units)) > 1) stop("unsupported data structure: the sampling units of all datasets must be the same.")
      if (length(unique(sampling_units)) == 1) n = unique(sampling_units)
      
      gamma = Reduce('+', data)
      gamma[gamma>1] = 1
      data_gamma_raw = gamma
      data_gamma_freq = c(n, rowSums(gamma))
      
      data_alpha_freq = sapply(data, rowSums) %>% c(n, .)
      
      
      data_2D = apply(sapply(data, rowSums), 2, function(x) c(n, x)) %>% as.data.frame
      
      ref_gamma = n
      ref_alpha = n
      
    }
    
    
    
    if (diversity == 'TD') {
      
      if (datatype == 'abundance') {
        
        gamma = estimate3D(as.numeric(data_gamma), diversity = 'TD', q = q, datatype = "abundance", base = "size", level = level, nboot = nboot) %>% arrange(., SC, Order.q)
        
        alpha = estimate3D(as.numeric(data_alpha), diversity = 'TD', q = q, datatype = "abundance", base = "size", level = level, nboot = nboot) %>% arrange(., SC, Order.q)
        
      }
      
      if (datatype == 'incidence_raw') {
        
        gamma = estimate3D(as.numeric(data_gamma_freq), diversity = 'TD', q = q, datatype = "incidence_freq", base = "size", level = level, nboot = nboot) %>% arrange(., SC, Order.q)
        
        alpha = estimate3D(data_alpha_freq, diversity = 'TD', q = q, datatype = "incidence_freq", base = "size", level = level, nboot = nboot) %>% arrange(., SC, Order.q)
        
      }
      
      se = cbind(gamma$s.e., alpha$s.e. / N)
      colnames(se) = c("gamma", "alpha")
      se = as.data.frame(se)
      # se[is.na(se)] = 0
      
      gamma = (cbind(Size = rep(level, each=length(q)), gamma[,-c(1,3,8,9)]) %>% 
                 mutate(Method = ifelse(Size>=ref_gamma, ifelse(Size == ref_gamma, 'Observed', 'Extrapolation'), 'Rarefaction'))
      )[,c(5,2,3,4,1)] %>% set_colnames(c('Estimate', 'Order.q', 'Method', 'SC', 'Size'))
      
      
      alpha = (cbind(Size = rep(level, each = length(q)), alpha[,-c(1,3,8,9)]) %>% 
                 mutate(Method = ifelse(Size >= ref_alpha, ifelse(Size == ref_alpha, 'Observed', 'Extrapolation'), 'Rarefaction'))
      )[,c(5,2,3,4,1)] %>% set_colnames(c('Estimate', 'Order.q', 'Method', 'SC', 'Size'))
      
      alpha$Estimate = alpha$Estimate / N
      
    }
    
    if (diversity == 'PD') {
      
      if (datatype == 'abundance') {
        
        aL = iNEXT.3D:::phyBranchAL_Abu(phylo = PDtree, data = data_gamma, rootExtend = T, refT = reft)
        aL$treeNabu$branch.length = aL$BLbyT[,1]
        aL_table_gamma = aL$treeNabu %>% select(branch.abun, branch.length, tgroup)
        
        gamma = iNEXT.3D:::PhD.m.est(ai = aL_table_gamma$branch.abun, Lis = as.matrix(aL_table_gamma$branch.length), m = level, q = q, nt = n, cal = "PD") %>% t %>% as.data.frame %>% 
          set_colnames(q) %>% gather(Order.q, Estimate) %>% 
          mutate(Coverage_real = rep(iNEXT.3D:::Coverage(data_gamma, "abundance", level), length(q)), Size = rep(level, length(q)), Size = rep(level, length(q)))
        
        
        aL_table_alpha = c()
        
        for (i in 1:N){
          
          x = data[data[,i]>0,i]
          names(x) = rownames(data)[data[,i]>0]
          
          aL = iNEXT.3D:::phyBranchAL_Abu(phylo = PDtree, data = x, rootExtend = T, refT = reft)
          aL$treeNabu$branch.length = aL$BLbyT[,1]
          aL_table = aL$treeNabu %>% select(branch.abun, branch.length, tgroup)
          
          aL_table_alpha = rbind(aL_table_alpha, aL_table)
          
        }
        
        
        qPDm = iNEXT.3D:::PhD.m.est(ai = aL_table_alpha$branch.abun, Lis = as.matrix(aL_table_alpha$branch.length), m = level, q = q, nt = n, cal = "PD")
        qPDm = qPDm/N
        alpha = qPDm %>% t %>% as.data.frame %>% 
          set_colnames(q) %>% gather(Order.q, Estimate) %>% 
          mutate(Coverage_real = rep(iNEXT.3D:::Coverage(data_alpha, "abundance", level), length(q)), Size = rep(level, length(q)))
        
      }
      
      if (datatype == 'incidence_raw') {
        
        aL = iNEXT.3D:::phyBranchAL_Inc(phylo = PDtree, data = as.matrix(data_gamma_raw), datatype = "incidence_raw", refT = reft, rootExtend = T)
        aL$treeNabu$branch.length = aL$BLbyT[,1]
        aL_table_gamma = aL$treeNabu %>% select(branch.abun, branch.length, tgroup)
        gamma = iNEXT.3D:::PhD.m.est(ai = aL_table_gamma$branch.abun, Lis=as.matrix(aL_table_gamma$branch.length), m = level, q = q, nt = n, cal = "PD") %>% t %>% as.data.frame %>% 
          set_colnames(q) %>% gather(Order.q, Estimate) %>% 
          mutate(Coverage_real = rep(iNEXT.3D:::Coverage(data_gamma_freq, "incidence_freq", level), length(q)), Size = rep(level, length(q)))
        
        aL_table_alpha = c()
        
        for (i in 1:N){
          
          x = data[[i]]
          
          aL = iNEXT.3D:::phyBranchAL_Inc(phylo = PDtree, data = x, datatype = "incidence_raw", rootExtend = T, refT = reft)
          aL$treeNabu$branch.length = aL$BLbyT[,1]
          aL_table = aL$treeNabu %>% select(branch.abun, branch.length, tgroup)
          
          aL_table_alpha = rbind(aL_table_alpha, aL_table)
          
        }
        
        alpha = (iNEXT.3D:::PhD.m.est(ai = aL_table_alpha$branch.abun, Lis = as.matrix(aL_table_alpha$branch.length), m = level, q = q, nt = n, cal = "PD")/N) %>% t %>% as.data.frame %>% 
          set_colnames(q) %>% gather(Order.q, Estimate) %>% 
          mutate(Coverage_real = rep(iNEXT.3D:::Coverage(data_alpha_freq, "incidence_freq", level), length(q)), Size = rep(level, length(q)))
        
        
      }
      
      gamma = (gamma %>% 
                 mutate(Method = ifelse(Size >= ref_gamma, ifelse(Size == ref_gamma, 'Observed', 'Extrapolation'), 'Rarefaction')))[,c(2,1,5,3,4)] %>% 
        set_colnames(c('Estimate', 'Order.q', 'Method', 'SC', 'Size'))
      
      gamma$Order.q = as.numeric(gamma$Order.q)
      
      
      alpha = (alpha %>% 
                 mutate(Method = ifelse(Size >= ref_alpha, ifelse(Size == ref_alpha, 'Observed', 'Extrapolation'), 'Rarefaction')))[,c(2,1,5,3,4)] %>% 
        set_colnames(c('Estimate', 'Order.q', 'Method', 'SC', 'Size'))
      
      alpha$Order.q = as.numeric(alpha$Order.q)
      
      if (PDtype == 'meanPD') {
        gamma$Estimate = gamma$Estimate/reft
        alpha$Estimate = alpha$Estimate/reft
      }
      
      # beta = alpha
      # beta$Estimate = gamma$Estimate/alpha$Estimate
      # 
      # C = beta %>% mutate(Estimate = ifelse(Order.q == 1, log(Estimate)/log(N), (Estimate^(1 - Order.q) - 1)/(N^(1 - Order.q) - 1)))
      # U = beta %>% mutate(Estimate = ifelse(Order.q == 1, log(Estimate)/log(N), (Estimate^(Order.q - 1) - 1)/(N^(Order.q - 1) - 1)))
      # V = beta %>% mutate(Estimate = (Estimate - 1)/(N - 1))
      # S = beta %>% mutate(Estimate = (1/Estimate - 1)/(1/N - 1))
      
      if(nboot>1){
        
        # cl = makeCluster(cluster_numbers)
        # clusterExport(cl, c("bootstrap_population_multiple_assemblage","data","data_gamma", 'data_gamma_freq',"level","N",'under_max_alpha',
        #                     'datatype', 'data_2D'))
        # clusterEvalQ(cl, library(tidyverse, magrittr))
        
        # plan(sequential)
        # plan(multiprocess)
        
        # se = parSapply(cl, 1:nboot, function(i){
        
        # start = Sys.time()
        se = future_lapply(1:nboot, function(i){
          
          if (datatype == 'abundance') {
            
            tree_bt = PDtree
            
            bootstrap_population = bootstrap_population_multiple_assemblage(data, data_gamma, 'abundance')
            p_bt = bootstrap_population
            unseen_p = p_bt[-(1:nrow(data)),] %>% matrix(ncol = ncol(data))
            
            if ( nrow(p_bt) > nrow(data) & sum(unseen_p) > 0 ){
              
              unseen = unseen_p[which(rowSums(unseen_p) > 0),]
              unseen = matrix(unseen, ncol = ncol(unseen_p), byrow = T)
              p_bt = rbind(p_bt[(1:nrow(data)),], unseen)
              unseen_name = sapply(1:nrow(unseen), function(i) paste0('unseen_', i))
              rownames(p_bt) = c(rownames(data), unseen_name)
              
              bootstrap_sample = sapply(1:ncol(data), function(k) rmultinom(n = 1, size = sum(data[,k]), prob = p_bt[,k]))
              x_bt = bootstrap_sample
              
              rownames(x_bt) = rownames(p_bt)
              
              if ( sum(x_bt[-(1:nrow(data)),])>0 ){
                
                g0_hat = apply(data, 2, function(x){
                  
                  n = sum(x)
                  f1 = sum(x == 1)
                  f2 = sum(x == 2)
                  
                  aL = iNEXT.3D:::phyBranchAL_Abu(phylo = PDtree, data = x, rootExtend = T, refT = reft)
                  
                  aL$treeNabu$branch.length = aL$BLbyT[,1]
                  aL = aL$treeNabu %>% select(branch.abun,branch.length)
                  g1 = aL$branch.length[aL$branch.abun == 1] %>% sum
                  g2 = aL$branch.length[aL$branch.abun == 2] %>% sum
                  g0_hat = ifelse( g2 > ((g1*f2)/(2*f1)) , ((n-1)/n)*(g1^2/(2*g2)) , ((n-1)/n)*(g1*(f1-1)/(2*(f2+1))) )
                  g0_hat
                  
                })
                
                te = (x_bt[1:nrow(data),]*(data == 0))>0
                used_length = sapply(1:ncol(data), function(i) { 
                  
                  if (sum(te[,i]) == 0) return(0) else {
                    
                    iNEXT.3D:::phyBranchAL_Abu(phylo = PDtree, data = x_bt[1:nrow(data),i], rootExtend = T, refT = reft)$treeNabu %>%
                      subset(label %in% names(which(te[,i] == TRUE))) %>% select(branch.length) %>% sum
                    
                  }
                  
                })
                
                g0_hat = g0_hat - used_length
                g0_hat[g0_hat < 0] = 0
                
                unseen_sample = x_bt[-(1:nrow(data)),]
                if (is.vector(unseen_sample)) unseen_sample = matrix(unseen_sample, ncol = ncol(x_bt), byrow = T)
                
                L0_hat = sapply(1:length(g0_hat), function(i) if(sum(unseen_sample[,i] > 0) > 0) (g0_hat[i] / nrow(unseen)) else 0 )
                
                L0_hat = rowSums((matrix(L0_hat, nrow(unseen_sample), ncol(unseen_sample), byrow = T) * unseen_sample)) / rowSums(unseen_sample)
                L0_hat[which(rowSums(unseen_sample) == 0)] = 0
                
                for (i in 1:length(L0_hat)){
                  
                  tip = list(edge = matrix(c(2,1),1,2),
                             tip.label = unseen_name[i],
                             edge.length = L0_hat[i],
                             Nnode = 1)
                  class(tip) = "phylo"
                  
                  tree_bt = tree_bt + tip
                  
                }
                
              } else {
                
                x_bt = x_bt[1:nrow(data),]
                p_bt = p_bt[1:nrow(data),]
                
              }
              
            } else {
              
              p_bt = p_bt[1:nrow(data),]
              x_bt = sapply(1:ncol(data), function(k) rmultinom(n = 1, size = sum(data[,k]), prob = p_bt[,k]))
              rownames(x_bt) = rownames(data)
              
            }
            
            bootstrap_data_gamma = rowSums(x_bt)
            bootstrap_data_gamma = bootstrap_data_gamma[bootstrap_data_gamma>0]
            bootstrap_data_alpha = as.matrix(x_bt) %>% as.vector
            bootstrap_data_alpha = bootstrap_data_alpha[bootstrap_data_alpha>0]
            
            aL = iNEXT.3D:::phyBranchAL_Abu(phylo = tree_bt, data = bootstrap_data_gamma, rootExtend = T, refT = reft)
            aL$treeNabu$branch.length = aL$BLbyT[,1]
            aL_table_gamma = aL$treeNabu %>% select(branch.abun, branch.length, tgroup)
            
            gamma = as.vector(iNEXT.3D:::PhD.m.est(ai = aL_table_gamma$branch.abun, Lis = as.matrix(aL_table_gamma$branch.length), m = level, q = q, nt = n, cal = "PD") %>% t)
            
            
            aL_table_alpha = c()
            
            for (i in 1:N){
              
              # x = x_bt[x_bt[,i]>0,i]
              # names(x) = rownames(p_bt)[x_bt[,i]>0]
              
              x = x_bt[,i]
              names(x) = rownames(p_bt)
              x = x[x_bt[,i]>0]
              
              aL = iNEXT.3D:::phyBranchAL_Abu(phylo = tree_bt, data = x, rootExtend = T, refT = reft)
              aL$treeNabu$branch.length = aL$BLbyT[,1]
              aL_table = aL$treeNabu %>% select(branch.abun, branch.length, tgroup)
              
              aL_table_alpha = rbind(aL_table_alpha, aL_table)
              
            }
            
            alpha = as.vector((iNEXT.3D:::PhD.m.est(ai = aL_table_alpha$branch.abun, Lis = as.matrix(aL_table_alpha$branch.length), m = level, q = q, nt = n, cal = "PD")/N) %>% t)
            
          }
          
          if (datatype == 'incidence_raw') {
            
            tree_bt = PDtree
            
            bootstrap_population = bootstrap_population_multiple_assemblage(data_2D, data_gamma_freq, 'incidence')
            p_bt = bootstrap_population
            unseen_p = p_bt[-(1:nrow(data[[1]])),] %>% matrix(ncol=N)
            
            if ( nrow(p_bt) > nrow(data[[1]]) & sum(unseen_p)>0 ){
              
              unseen = unseen_p[which(rowSums(unseen_p)>0),]
              unseen = matrix(unseen, ncol = ncol(unseen_p), byrow = T)
              p_bt = rbind(p_bt[(1:nrow(data[[1]])),], unseen)
              unseen_name = sapply(1:nrow(unseen), function(i) paste0('unseen_', i))
              rownames(p_bt) = c(rownames(data[[1]]), unseen_name)
              
              raw = lapply(1:ncol(p_bt), function(j){
                
                lapply(1:nrow(p_bt), function(i) rbinom(n=n, size=1, prob=p_bt[i,j])) %>% do.call(rbind,.)
                
              })
              
              for (i in 1:length(raw)) rownames(raw[[i]]) = rownames(p_bt)
              
              if ( lapply(1:length(raw), function(i) raw[[i]][-(1:nrow(data[[1]])),]) %>% do.call(sum,.)>0 ){
                
                R0_hat = sapply(data, function(x){
                  
                  nT = ncol(x)
                  Q1 = sum(rowSums(x)==1)
                  Q2 = sum(rowSums(x)==2)
                  
                  aL = iNEXT.3D:::phyBranchAL_Inc(phylo = PDtree, data = x, datatype = "incidence_raw", rootExtend = T, refT = reft)
                  
                  aL$treeNabu$branch.length = aL$BLbyT[,1]
                  aL = aL$treeNabu %>% select(branch.abun,branch.length)
                  R1 = aL$branch.length[aL$branch.abun == 1] %>% sum
                  R2 = aL$branch.length[aL$branch.abun == 2] %>% sum
                  R0_hat = ifelse( R2>((R1*Q2)/(2*Q1)) , ((nT-1)/nT)*(R1^2/(2*R2)) , ((nT-1)/nT)*(R1*(Q1-1)/(2*(Q2+1))) )
                  R0_hat
                  
                })
                
                te = (sapply(raw, rowSums)[1:nrow(data[[1]]),]*(sapply(data, rowSums) == 0)) > 0
                used_length = sapply(1:N, function(i) {
                  
                  if (sum(te[,i]) == 0) return(0) else {
                    
                    iNEXT.3D:::phyBranchAL_Inc(phylo = PDtree, data = raw[[i]][1:nrow(data[[1]]),], datatype = "incidence_raw", rootExtend = T, refT = reft)$treeNabu %>%
                      subset(label %in% names(which(te[,i] == TRUE))) %>% select(branch.length) %>% sum
                    
                  }
                  
                })
                
                R0_hat = R0_hat - used_length
                R0_hat[R0_hat < 0] = 0
                
                unseen_sample = sapply(raw, rowSums)[-(1:nrow(data[[1]])),]
                if (is.vector(unseen_sample)) unseen_sample = matrix(unseen_sample, ncol = N, byrow = T)
                
                L0_hat = sapply(1:length(R0_hat), function(i) if(sum(unseen_sample[,i] > 0) > 0) (R0_hat[i] / nrow(unseen)) else 0 )
                
                L0_hat = rowSums((matrix(L0_hat, nrow(unseen_sample), ncol(unseen_sample), byrow = T) * unseen_sample)) / rowSums(unseen_sample)
                L0_hat[which(rowSums(unseen_sample) == 0)] = 0
                
                for (i in 1:length(L0_hat)){
                  
                  tip = list(edge = matrix(c(2,1), 1, 2),
                             tip.label = unseen_name[i],
                             edge.length = L0_hat[i],
                             Nnode = 1)
                  class(tip) = "phylo"
                  
                  tree_bt = tree_bt + tip
                  
                }
                
              } else raw = lapply(raw, function(i) i[1:nrow(data[[1]]),])
              
            } else {
              
              p_bt = p_bt[1:nrow(data[[1]]),]
              raw = lapply(1:ncol(p_bt), function(j){
                
                lapply(1:nrow(p_bt), function(i) rbinom(n = n, size = 1, prob = p_bt[i,j])) %>% do.call(rbind,.)
                
              })
              
              for (i in 1:length(raw)) rownames(raw[[i]]) = rownames(p_bt)
              
            }
            
            gamma = Reduce('+', raw)
            gamma[gamma > 1] = 1
            bootstrap_data_gamma_raw = gamma
            bootstrap_data_gamma_freq = c(n, rowSums(gamma))
            
            bootstrap_data_alpha_freq = sapply(raw, rowSums) %>% c(n, .)
            
            bootstrap_data_gamma_freq = bootstrap_data_gamma_freq[bootstrap_data_gamma_freq > 0]
            bootstrap_data_alpha_freq = bootstrap_data_alpha_freq[bootstrap_data_alpha_freq > 0]
            
            aL = iNEXT.3D:::phyBranchAL_Inc(phylo = tree_bt, data = bootstrap_data_gamma_raw, datatype = "incidence_raw", rootExtend = T, refT = reft)
            aL$treeNabu$branch.length = aL$BLbyT[,1]
            aL_table_gamma = aL$treeNabu %>% select(branch.abun, branch.length, tgroup)
            
            gamma = as.vector(iNEXT.3D:::PhD.m.est(ai = aL_table_gamma$branch.abun, Lis = as.matrix(aL_table_gamma$branch.length), m = level, q = q, nt = n, cal="PD") %>% t)
            
            
            aL_table_alpha = c()
            
            for (i in 1:N){
              
              x = raw[[i]]
              
              aL = iNEXT.3D:::phyBranchAL_Inc(phylo = tree_bt, data = x, datatype = "incidence_raw", rootExtend = T, refT = reft)
              aL$treeNabu$branch.length = aL$BLbyT[,1]
              aL_table = aL$treeNabu %>% select(branch.abun, branch.length, tgroup)
              
              aL_table_alpha = rbind(aL_table_alpha, aL_table)
              
            }
            
            alpha = as.vector((iNEXT.3D:::PhD.m.est(ai = aL_table_alpha$branch.abun, Lis = as.matrix(aL_table_alpha$branch.length), m = level, q=q, nt = n, cal = "PD")/N) %>% t)
            
          }
          
          if (PDtype == 'meanPD') {
            gamma = gamma/reft
            alpha = alpha/reft
          } 
          
          # beta = gamma/alpha
          # 
          # Order.q = rep(q, each = length(level))
          # 
          # beta = data.frame(Estimate = beta, Order.q)
          # 
          # C = (beta %>% mutate(Estimate = ifelse(Order.q == 1, log(Estimate)/log(N), (Estimate^(1 - Order.q) - 1)/(N^(1 - Order.q) - 1))))$Estimate
          # U = (beta %>% mutate(Estimate = ifelse(Order.q == 1, log(Estimate)/log(N), (Estimate^(Order.q - 1) - 1)/(N^(Order.q - 1) - 1))))$Estimate
          # V = (beta %>% mutate(Estimate = (Estimate - 1)/(N - 1)))$Estimate
          # S = (beta %>% mutate(Estimate = (1/Estimate - 1)/(1/N - 1)))$Estimate
          # 
          # beta = beta$Estimate
          # 
          # cbind(gamma, alpha, beta, C, U, V, S) %>% as.matrix
          cbind(gamma, alpha) %>% as.matrix
          
          # }, simplify = "array") %>% apply(., 1:2, sd) %>% data.frame
        }, future.seed = TRUE) %>% abind(along = 3) %>% apply(1:2, sd)
        # end = Sys.time()
        # end - start
        
        # stopCluster(cl) 
        # plan(sequential)
        
      } else {
        
        # se = matrix(0, ncol = 7, nrow = nrow(beta))
        # colnames(se) = c("gamma", "alpha", "beta", "C", "U", 'V', 'S')
        # se = as.data.frame(se)
        
        se = matrix(NA, ncol = 2, nrow = nrow(gamma))
        colnames(se) = c("gamma", "alpha")
        se = as.data.frame(se)
        
      }
      
    }
    
    if (diversity == 'FD') {
      
      FDdistM = as.matrix(FDdistM)
      
      FD_by_tau = function(data, distM, tau, datatype, m_gamma, m_alpha) {
        
        if (datatype == 'abundance') {
          
          zik = data
          zik = zik[rowSums(data)>0,]
          
          dij = distM
          dij = dij[rowSums(data)>0, rowSums(data)>0]
          
          if (tau == 0) {
            dij[dij > 0] <- 1
            aik = (1 - dij/1) %*% as.matrix(zik)
          } else {
            dij[which(dij>tau, arr.ind = T)] = tau
            aik = (1 - dij/tau) %*% as.matrix(zik)
          }
          positive_id = rowSums(aik)>0
          
          
          gamma_x = rowSums(zik)[positive_id]
          gamma_a = rowSums(aik)[positive_id]
          gamma_v = gamma_x/gamma_a
          # gamma_a = ifelse(gamma_a < 1, 1, round(gamma_a))
          
          ai_vi_gamma = list(ai = data.frame(gamma_a), vi = data.frame(gamma_v))
          
          gamma = FD.m.est_0(ai_vi_gamma, m_gamma, q, n) %>% as.vector
          # gamma = iNEXT.3D:::FD.m.est(ai_vi_gamma, m_gamma, q, n) %>% as.vector
          
          
          alpha_x = as.vector(as.matrix(zik))
          alpha_a = as.vector(aik)
          # alpha_a = ifelse(alpha_a < 1, 1, round(alpha_a))
          
          alpha_v = alpha_x/alpha_a
          alpha_v = rep(gamma_v,N)
          
          alpha_v = alpha_v[alpha_a>0]
          alpha_a = alpha_a[alpha_a>0]
          
          ai_vi_alpha = list(ai = data.frame(alpha_a), vi = data.frame(alpha_v))
          
          alpha = (FD.m.est_0(ai_vi_alpha, m_gamma, q, n)/N) %>% as.vector
          
        }
        
        if (datatype == 'incidence_raw') {
          
          data_gamma_freq = data$data_gamma_freq
          data_2D = data$data_2D
          
          gamma_Y = data_gamma_freq[-1]
          
          dij = distM
          dij = dij[gamma_Y > 0, gamma_Y > 0]
          gamma_Y = gamma_Y[gamma_Y > 0]
          
          if (tau == 0) {
            dij[dij > 0] <- 1
            gamma_a = (1 - dij/1) %*% as.matrix(gamma_Y)
          } else {
            dij[which(dij>tau, arr.ind = T)] = tau
            gamma_a = (1 - dij/tau) %*% as.matrix(gamma_Y)
          }
          
          gamma_a[gamma_a > n] = n
          gamma_v = gamma_Y/gamma_a
          
          # gamma_a = ifelse(gamma_a < 1, 1, round(gamma_a))
          # gamma_a[gamma_a > n] = n
          
          ai_vi_gamma = list(ai = data.frame(gamma_a), vi = data.frame(gamma_v))
          
          gamma = FD.m.est_0(ai_vi_gamma, m_gamma, q, n) %>% as.vector
          # gamma = iNEXT.3D:::FD.m.est(ai_vi_gamma, m_gamma, q, n) %>% as.vector
          
          
          alpha_Y = data_2D[-1,]
          
          dij = distM
          dij = dij[rowSums(data_2D[-1,]) > 0, rowSums(data_2D[-1,])>0]
          alpha_Y = alpha_Y[rowSums(data_2D[-1,])>0,]
          
          
          if (tau == 0) {
            dij[dij > 0] <- 1
            alpha_a = (1 - dij/1) %*% as.matrix(alpha_Y)
          } else {
            dij[which(dij>tau, arr.ind = T)] = tau
            alpha_a = (1 - dij/tau) %*% as.matrix(alpha_Y)
          }
          
          # alpha_a = ifelse(alpha_a < 1, 1, round(alpha_a))
          alpha_a[alpha_a > n] = n
          alpha_a = as.vector(alpha_a)
          
          alpha_v = rep(gamma_v, N)
          alpha_v = alpha_v[alpha_a > 0]
          alpha_a = alpha_a[alpha_a > 0]
          
          ai_vi_alpha = list(ai = data.frame(alpha_a), vi = data.frame(alpha_v))
          alpha = (FD.m.est_0(ai_vi_alpha, m_gamma, q, n)/N) %>% as.vector
          
        }
        
        return(data.frame(gamma,alpha))
        
      }
      
      if (FDtype == 'tau_value'){
        
        if (datatype == 'abundance') {
          
          FDdistM = FDdistM[rownames(FDdistM) %in% rownames(data), colnames(FDdistM) %in% rownames(data)]
          order_sp <- match(rownames(data),rownames(FDdistM))
          FDdistM <- FDdistM[order_sp,order_sp]
          
          output = FD_by_tau(data, FDdistM, FDtau, datatype = 'abundance', m_gamma = level, m_alpha = level)
          gamma = output$gamma
          alpha = output$alpha
          
          gamma = data.frame(Size = level, gamma) %>% 
            mutate(Method = ifelse(Size >= ref_gamma, ifelse(Size == ref_gamma, 'Observed', 'Extrapolation'), 'Rarefaction'),
                   Order.q = rep(q, each = length(level)), Coverage_real = rep(iNEXT.3D:::Coverage(data_gamma, "abundance", level), length(q)))
          
          alpha = data.frame(Size = level, alpha) %>% 
            mutate(Method = ifelse(Size>=ref_alpha, ifelse(Size == ref_alpha, 'Observed', 'Extrapolation'), 'Rarefaction'),
                   Order.q = rep(q, each = length(level)), Coverage_real = rep(iNEXT.3D:::Coverage(data_alpha, "abundance", level), length(q)))
          
          beta = alpha
          beta$alpha = gamma$gamma/alpha$alpha
          
        }
        
        if (datatype == 'incidence_raw') {
          
          FDdistM = FDdistM[rownames(FDdistM) %in% names(data_gamma_freq)[-1], colnames(FDdistM) %in% names(data_gamma_freq)[-1]]
          order_sp <- match(names(data_gamma_freq)[-1],rownames(FDdistM))
          FDdistM <- FDdistM[order_sp,order_sp]
          
          output = FD_by_tau(list(data_gamma_freq = data_gamma_freq, data_2D = data_2D), FDdistM, FDtau, datatype='incidence_raw', m_gamma=level, m_alpha=level)
          gamma = output$gamma
          alpha = output$alpha
          
          gamma = data.frame(Size = level, gamma) %>% 
            mutate(Method = ifelse(Size >= ref_gamma, ifelse(Size == ref_gamma, 'Observed', 'Extrapolation'), 'Rarefaction'),
                   Order.q = rep(q, each = length(level)), Coverage_real = rep(iNEXT.3D:::Coverage(data_gamma_freq, "incidence_freq", level), length(q)))
          
          alpha = data.frame(Size = level, alpha) %>% 
            mutate(Method = ifelse(Size>=ref_alpha, ifelse(Size == ref_alpha, 'Observed', 'Extrapolation'), 'Rarefaction'),
                   Order.q = rep(q, each = length(level)), Coverage_real = rep(iNEXT.3D:::Coverage(data_alpha_freq, "incidence_freq", level), length(q)))
          
          beta = alpha
          beta$alpha = gamma$gamma/alpha$alpha
          
        }
        
      }
      
      if (FDtype == 'AUC'){
        
        cut = seq(0.00000001, 1, length.out = FDcut_number)
        width = diff(cut)
        
        if (datatype == 'abundance') {
          
          FDdistM = FDdistM[rownames(FDdistM) %in% rownames(data), colnames(FDdistM) %in% rownames(data)]
          order_sp <- match(rownames(data),rownames(FDdistM))
          FDdistM <- FDdistM[order_sp,order_sp]
          
          gamma_alpha_over_tau = lapply(cut, function(tau) {
            
            FD_by_tau(data, FDdistM, tau, datatype = 'abundance', m_gamma = level, m_alpha = level)
            
          })
          
          gamma_over_tau = sapply(gamma_alpha_over_tau, function(x) x$gamma)
          
          left_limit  = apply(gamma_over_tau, 1, function(x) x[-FDcut_number]*width)
          right_limit = apply(gamma_over_tau, 1, function(x) x[-1]*width)
          
          gamma = colSums((left_limit + right_limit)/2)
          
          alpha_over_tau = sapply(gamma_alpha_over_tau, function(x) x$alpha)
          
          left_limit  = apply(alpha_over_tau, 1, function(x) x[-FDcut_number]*width)
          right_limit = apply(alpha_over_tau, 1, function(x) x[-1]*width)
          
          alpha = colSums((left_limit + right_limit)/2)
          
          beta_over_tau = gamma_over_tau/alpha_over_tau
          
          left_limit  = apply(beta_over_tau, 1, function(x) x[-FDcut_number]*width)
          right_limit = apply(beta_over_tau, 1, function(x) x[-1]*width)
          
          beta = colSums((left_limit + right_limit)/2)
          
          gamma = data.frame(Size = level, gamma) %>% 
            mutate(Method = ifelse(Size >= ref_gamma, ifelse(Size == ref_gamma, 'Observed', 'Extrapolation'), 'Rarefaction'),
                   Order.q = rep(q, each = length(level)), Coverage_real = rep(iNEXT.3D:::Coverage(data_gamma, "abundance", level), length(q)))
          
          alpha = data.frame(Size = level, alpha) %>% 
            mutate(Method = ifelse(Size >= ref_alpha, ifelse(Size == ref_alpha, 'Observed', 'Extrapolation'), 'Rarefaction'),
                   Order.q = rep(q, each = length(level)), Coverage_real = rep(iNEXT.3D:::Coverage(data_alpha, "abundance", level), length(q)))
          
          # beta = data.frame(Size = level, beta) %>% 
          #   mutate(Method = ifelse(Size >= ref_alpha, ifelse(Size == ref_alpha, 'Observed', 'Extrapolation'), 'Rarefaction'),
          #          Order.q = rep(q, each = length(level)), Coverage_real = rep(iNEXT.3D:::Coverage(data_alpha, "abundance", level), length(q)))
          
        }
        
        if (datatype == 'incidence_raw') {
          
          FDdistM = FDdistM[rownames(FDdistM) %in% names(data_gamma_freq)[-1], colnames(FDdistM) %in% names(data_gamma_freq)[-1]]
          order_sp <- match(names(data_gamma_freq)[-1],rownames(FDdistM))
          FDdistM <- FDdistM[order_sp,order_sp]
          
          gamma_alpha_over_tau = lapply(cut, function(tau) {
            
            FD_by_tau(list(data_gamma_freq = data_gamma_freq, data_2D = data_2D), FDdistM, tau, datatype = 'incidence_raw', m_gamma = level, m_alpha = level)
            
          })
          
          gamma_over_tau = sapply(gamma_alpha_over_tau, function(x) x$gamma)
          
          left_limit  = apply(gamma_over_tau, 1, function(x) x[-FDcut_number]*width)
          right_limit = apply(gamma_over_tau, 1, function(x) x[-1]*width)
          
          gamma = colSums((left_limit + right_limit)/2)
          
          alpha_over_tau = sapply(gamma_alpha_over_tau, function(x) x$alpha)
          
          left_limit  = apply(alpha_over_tau, 1, function(x) x[-FDcut_number]*width)
          right_limit = apply(alpha_over_tau, 1, function(x) x[-1]*width)
          
          alpha = colSums((left_limit + right_limit)/2)
          
          beta_over_tau = gamma_over_tau/alpha_over_tau
          
          left_limit  = apply(beta_over_tau, 1, function(x) x[-FDcut_number]*width)
          right_limit = apply(beta_over_tau, 1, function(x) x[-1]*width)
          
          beta = colSums((left_limit + right_limit)/2)
          
          gamma = data.frame(Size = level, gamma) %>% 
            mutate(Method = ifelse(Size >= ref_gamma, ifelse(Size == ref_gamma, 'Observed', 'Extrapolation'), 'Rarefaction'),
                   Order.q = rep(q, each = length(level)), Coverage_real = rep(iNEXT.3D:::Coverage(data_gamma_freq, "incidence_freq", level), length(q)))
          
          alpha = data.frame(Size = level, alpha) %>% 
            mutate(Method = ifelse(Size >= ref_alpha, ifelse(Size == ref_alpha, 'Observed', 'Extrapolation'), 'Rarefaction'),
                   Order.q = rep(q, each = length(level)), Coverage_real = rep(iNEXT.3D:::Coverage(data_alpha_freq, "incidence_freq", level), length(q)))
          
          # beta = data.frame(Size = level, beta) %>% 
          #   mutate(Method = ifelse(Size >= ref_alpha, ifelse(Size == ref_alpha, 'Observed', 'Extrapolation'), 'Rarefaction'),
          #          Order.q = rep(q, each = length(level)), Coverage_real = rep(iNEXT.3D:::Coverage(data_alpha_freq, "incidence_freq", level), length(q)))
          
        }
        
      }
      
      gamma = gamma[,c(2,4,3,5,1)] %>% set_colnames(c('Estimate', 'Order.q', 'Method', 'SC', 'Size'))
      
      
      alpha = alpha[,c(2,4,3,5,1)] %>% set_colnames(c('Estimate', 'Order.q', 'Method', 'SC', 'Size'))
      
      
      # beta = beta[,c(2,4,3,5,1)] %>% set_colnames(c('Estimate', 'Order.q', 'Method', 'SC', 'Size'))
      # 
      # C = beta %>% mutate(Estimate = ifelse(Order.q == 1, log(Estimate)/log(N), (Estimate^(1 - Order.q) - 1)/(N^(1 - Order.q) - 1)))
      # U = beta %>% mutate(Estimate = ifelse(Order.q == 1, log(Estimate)/log(N), (Estimate^(Order.q - 1) - 1)/(N^(Order.q - 1) - 1)))
      # V = beta %>% mutate(Estimate = (Estimate - 1)/(N - 1))
      # S = beta %>% mutate(Estimate = (1/Estimate - 1)/(1/N - 1))
      
      if(nboot > 1){
        
        # cl = makeCluster(cluster_numbers)
        # clusterExport(cl, c("bootstrap_population_multiple_assemblage","data","data_gamma", 'data_gamma_freq',"level","N",'under_max_alpha',
        #                     'datatype', 'data_2D'))
        # clusterEvalQ(cl, library(tidyverse, magrittr))
        
        # plan(sequential)
        # plan(multiprocess, workers=7)
        
        # se = parSapply(cl, 1:nboot, function(i){
        
        # start = Sys.time()
        se = future_lapply(1:nboot, function(i){
          
          if (datatype == 'abundance') {
            
            p_bt = bootstrap_population_multiple_assemblage(data, data_gamma, 'abundance')
            f0_hat = nrow(p_bt) - nrow(data)
            
            distance_matrix_bt = Bootstrap_distance_matrix(rowSums(data), FDdistM, f0_hat, 'abundance')
            
            data_bt = sapply(1:ncol(data), function(k) rmultinom(n = 1, size = sum(data[,k]), prob = p_bt[,k]))
            
            data_gamma = rowSums(data_bt)
            data_gamma = data_gamma[data_gamma>0]
            data_alpha = as.matrix(data_bt) %>% as.vector
            
            if (FDtype == 'tau_value'){
              
              output = FD_by_tau(data_bt, distance_matrix_bt, FDtau, datatype='abundance', m_gamma = level, m_alpha = level)
              gamma = output$gamma
              alpha = output$alpha
              
              # beta=gamma/alpha
              
            }
            
            if (FDtype == 'AUC'){
              
              gamma_alpha_over_tau = lapply(cut, function(tau) {
                
                FD_by_tau(data_bt, distance_matrix_bt, tau, datatype = 'abundance', m_gamma = level, m_alpha = level)
                
              })
              
              gamma_over_tau = sapply(gamma_alpha_over_tau, function(x) x$gamma)
              
              left_limit  = apply(gamma_over_tau, 1, function(x) x[-FDcut_number]*width)
              right_limit = apply(gamma_over_tau, 1, function(x) x[-1]*width)
              
              gamma = colSums((left_limit + right_limit)/2)
              
              alpha_over_tau = sapply(gamma_alpha_over_tau, function(x) x$alpha)
              
              left_limit  = apply(alpha_over_tau, 1, function(x) x[-FDcut_number]*width)
              right_limit = apply(alpha_over_tau, 1, function(x) x[-1]*width)
              
              alpha = colSums((left_limit + right_limit)/2)
              
              # beta_over_tau = gamma_over_tau/alpha_over_tau
              # 
              # left_limit  = apply(beta_over_tau, 1, function(x) x[-FDcut_number]*width)
              # right_limit = apply(beta_over_tau, 1, function(x) x[-1]*width)
              # 
              # beta = colSums((left_limit + right_limit)/2)
              
            }
            
          }
          
          if (datatype == 'incidence_raw') {
            
            p_bt = bootstrap_population_multiple_assemblage(data_2D, data_gamma_freq, 'incidence')
            f0_hat = nrow(p_bt) - nrow(data_2D[-1,])
            
            distance_matrix_bt = Bootstrap_distance_matrix(c(n,rowSums(data_gamma_raw)), FDdistM, f0_hat, 'incidence_freq')
            
            # p_bt = p_bt[rowSums(p_bt)>0,]
            
            raw = lapply(1:ncol(p_bt), function(j){
              
              lapply(1:nrow(p_bt), function(i) rbinom(n = n, size = 1, prob = p_bt[i,j])) %>% do.call(rbind,.)
              
            })
            
            gamma = Reduce('+', raw)
            gamma[gamma>1] = 1
            data_gamma_raw_bt = gamma
            data_gamma_freq_bt = c(n, rowSums(gamma))
            
            data_alpha_freq_bt = sapply(raw, rowSums) %>% c(n, .)
            
            # data_gamma_freq_bt = data_gamma_freq_bt[data_gamma_freq_bt > 0]
            # data_alpha_freq_bt = data_alpha_freq_bt[data_alpha_freq_bt > 0]
            
            data_2D_bt = apply(sapply(raw, rowSums), 2, function(x) c(n, x)) %>% as.data.frame
            
            if (FDtype == 'tau_value'){
              
              output = FD_by_tau(list(data_gamma_freq = data_gamma_freq_bt, data_2D = data_2D_bt), distance_matrix_bt, FDtau, datatype = 'incidence_raw', m_gamma = level, m_alpha = level)
              gamma = output$gamma
              alpha = output$alpha
              
              # beta = gamma/alpha
              
            }
            
            if (FDtype == 'AUC'){
              
              gamma_alpha_over_tau = lapply(cut, function(tau) {
                
                FD_by_tau(list(data_gamma_freq = data_gamma_freq_bt, data_2D = data_2D_bt), distance_matrix_bt, tau, datatype = 'incidence_raw', m_gamma = level, m_alpha = level)
                
              })
              
              gamma_over_tau = sapply(gamma_alpha_over_tau, function(x) x$gamma)
              
              left_limit  = apply(gamma_over_tau, 1, function(x) x[-FDcut_number]*width)
              right_limit = apply(gamma_over_tau, 1, function(x) x[-1]*width)
              
              gamma = colSums((left_limit + right_limit)/2)
              
              alpha_over_tau = sapply(gamma_alpha_over_tau, function(x) x$alpha)
              
              left_limit  = apply(alpha_over_tau, 1, function(x) x[-FDcut_number]*width)
              right_limit = apply(alpha_over_tau, 1, function(x) x[-1]*width)
              
              alpha = colSums((left_limit + right_limit)/2)
              
              # beta_over_tau = gamma_over_tau/alpha_over_tau
              # 
              # left_limit  = apply(beta_over_tau, 1, function(x) x[-FDcut_number]*width)
              # right_limit = apply(beta_over_tau, 1, function(x) x[-1]*width)
              # 
              # beta = colSums((left_limit + right_limit)/2)
              
            }
            
          }
          
          # beta = gamma/alpha
          # 
          # Order.q = rep(q, each=length(level))
          # 
          # beta = data.frame(Estimate = beta, Order.q)
          # 
          # C = (beta %>% mutate(Estimate = ifelse(Order.q == 1, log(Estimate)/log(N), (Estimate^(1 - Order.q) - 1)/(N^(1 - Order.q) - 1))))$Estimate
          # U = (beta %>% mutate(Estimate = ifelse(Order.q == 1, log(Estimate)/log(N), (Estimate^(Order.q - 1) - 1)/(N^(Order.q - 1) - 1))))$Estimate
          # V = (beta %>% mutate(Estimate = (Estimate - 1)/(N - 1)))$Estimate
          # S = (beta %>% mutate(Estimate = (1/Estimate - 1)/(1/N - 1)))$Estimate
          # 
          # beta = beta$Estimate
          # 
          # cbind(gamma, alpha, beta, C, U, V, S) %>% as.matrix
          cbind(gamma, alpha) %>% as.matrix
          
          # }, simplify = "array") %>% apply(., 1:2, sd) %>% data.frame
        }, future.seed = TRUE) %>% abind(along = 3) %>% apply(1:2, sd)
        # end = Sys.time()
        # end - start
        
        # stopCluster(cl)
        # plan(sequential)
        
      } else {
        
        # se = matrix(0, ncol = 7, nrow = nrow(beta))
        # colnames(se) = c("gamma", "alpha", "beta", "C", "U", 'V', 'S')
        # se = as.data.frame(se)
        
        se = matrix(NA, ncol = 2, nrow = nrow(gamma))
        colnames(se) = c("gamma", "alpha")
        se = as.data.frame(se)
        
      }
      
    }
    
    
    se = as.data.frame(se)
    
    if (diversity == "TD") index = "TD"
    if (diversity == "PD" & PDtype == "PD") index = "PD"
    if (diversity == "PD" & PDtype == "meanPD") index = "meanPD"
    if (diversity == "FD" & FDtype == "tau_value") index = "FD_tau"
    if (diversity == "FD" & FDtype == "AUC") index = "FD_AUC"
    
    gamma = gamma %>% mutate(s.e. = se$gamma,
                             LCL = Estimate - tmp * se$gamma,
                             UCL = Estimate + tmp * se$gamma,
                             Dataset = dataset_name,
                             Diversity = index) %>% 
      arrange(Order.q, Size) %>% .[,c(9, 2, 5, 4, 1, 3, 6, 7, 8, 10)] %>% rename("Gamma" = "Estimate")
    
    alpha = alpha %>% mutate(s.e. = se$alpha,
                             LCL = Estimate - tmp * se$alpha,
                             UCL = Estimate + tmp * se$alpha,
                             Dataset = dataset_name,
                             Diversity = index) %>% 
      arrange(Order.q, Size) %>% .[,c(9, 2, 5, 4, 1, 3, 6, 7, 8, 10)] %>% rename("Alpha" = "Estimate")
    
    if (datatype == 'incidence_raw') {
      
      gamma = gamma %>% rename("nT" = "Size")
      alpha = alpha %>% rename("nT" = "Size")
    }
    
    if (diversity == "FD" & FDtype == "tau_value") {
      
      gamma = gamma %>% mutate(Tau = FDtau)
      
      alpha = alpha %>% mutate(Tau = FDtau)
    }
    
    list(gamma = gamma, alpha = alpha)
    
  }

  if (base == 'coverage') output = lapply(1:length(data_list), function(i) for_each_dataset(data = data_list[[i]], dataset_name = dataset_names[i], N = Ns[i], level = level[[i]]))
  if (base == 'size') output = lapply(1:length(data_list), function(i) for_each_dataset.size(data = data_list[[i]], dataset_name = dataset_names[i], N = Ns[i], level = level[[i]]))
  names(output) = dataset_names
  
  class(output) <- c("iNEXTbeta3D")
  return(output)
  
}



#' ggplot2 extension for an iNEXTbeta3D object
#' 
#' \code{ggiNEXTbeta3D}: the \code{\link[ggplot2]{ggplot}} extension for \code{\link{iNEXTbeta3D}} 
#' object to plot sample-size- and coverage-based rarefaction/extrapolation curves.
#' 
#' @param output the output from the function iNEXTbeta3D
#' @param type (argument only for \code{base = "coverage"}),\cr
#' \code{type = 'B'} for plotting the rarefaction and extrapolation sampling curves for gamma, alpha, and beta diversity;  \cr
#' \code{type = 'D'} for plotting the rarefaction and extrapolation sampling curves for four dissimilarity indices.
#' Skip the argument for plotting size-based rarefaction and extrapolation sampling curves for gamma and alpha diversity.
#' 
#' @return a figure for gamma, alpha, and beta diversity, or a figure for four dissimilarity indices for \code{base = "coverage"}; 
#' or a figure for gamma and alpha diversity when \code{base = "size"}.\cr
#' 
#' @examples
#' ## Taxonomic diversity for abundance data
#' # Coverage-based rarefaction and extrapolation sampling curves 
#' data(Brazil_rainforests)
#' output1c = iNEXTbeta3D(data = Brazil_rainforests, diversity = 'TD', 
#'                        datatype = 'abundance', base = "coverage", nboot = 30, conf = 0.95)
#' 
#' ggiNEXTbeta3D(output1c, type = 'B')
#' ggiNEXTbeta3D(output1c, type = 'D')
#'  
#' 
#' # Size-based rarefaction and extrapolation sampling curves 
#' data(Brazil_rainforests)
#' output1s = iNEXTbeta3D(data = Brazil_rainforests, diversity = 'TD', 
#'                        datatype = 'abundance', base = "size", nboot = 30, conf = 0.95)
#' 
#' ggiNEXTbeta3D(output1s)
#' 
#' 
#' ## Phylogenetic diversity for abundance data
#' # Coverage-based rarefaction and extrapolation sampling curves 
#' data(Brazil_rainforests)
#' data(Brazil_tree)
#' output2c = iNEXTbeta3D(data = Brazil_rainforests, diversity = 'PD', 
#'                        datatype = 'abundance', base = "coverage", nboot = 10, conf = 0.95, 
#'                        PDtree = Brazil_tree, PDreftime = NULL, PDtype = 'PD')
#' 
#' ggiNEXTbeta3D(output2c, type = 'B')
#' ggiNEXTbeta3D(output2c, type = 'D')
#' 
#' 
#' # Size-based rarefaction and extrapolation sampling curves 
#' data(Brazil_rainforests)
#' data(Brazil_tree)
#' output2s = iNEXTbeta3D(data = Brazil_rainforests, diversity = 'PD', 
#'                        datatype = 'abundance', base = "size", nboot = 10, conf = 0.95, 
#'                        PDtree = Brazil_tree, PDreftime = NULL, PDtype = 'PD')
#' 
#' ggiNEXTbeta3D(output2s)
#' 
#' 
#' ## Functional diversity for abundance data under a specified threshold level
#' # Coverage-based rarefaction and extrapolation sampling curves 
#' data(Brazil_rainforests)
#' data(Brazil_distM)
#' output3c = iNEXTbeta3D(data = Brazil_rainforests, diversity = 'FD', 
#'                        datatype = 'abundance', base = "coverage", nboot = 30, conf = 0.95, 
#'                        FDdistM = Brazil_distM, FDtype = 'tau_value', FDtau = NULL)
#' 
#' ggiNEXTbeta3D(output3c, type = 'B')
#' ggiNEXTbeta3D(output3c, type = 'D')
#' 
#' 
#' # Size-based rarefaction and extrapolation sampling curves 
#' data(Brazil_rainforests)
#' data(Brazil_distM)
#' output3s = iNEXTbeta3D(data = Brazil_rainforests, diversity = 'FD', 
#'                        datatype = 'abundance', base = "size", nboot = 30, conf = 0.95, 
#'                        FDdistM = Brazil_distM, FDtype = 'tau_value', FDtau = NULL)
#' 
#' ggiNEXTbeta3D(output3s)
#' 
#' 
#' ## Functional diversity for abundance data when all threshold levels from 0 to 1 are considered
#' # Coverage-based rarefaction and extrapolation sampling curves 
#' data(Brazil_rainforests)
#' data(Brazil_distM)
#' output4c = iNEXTbeta3D(data = Brazil_rainforests, diversity = 'FD', 
#'                        datatype = 'abundance', base = "coverage", nboot = 0, conf = 0.95, 
#'                        FDdistM = Brazil_distM, FDtype = 'AUC', FDcut_number = 30)
#' 
#' ggiNEXTbeta3D(output4c, type = 'B')
#' ggiNEXTbeta3D(output4c, type = 'D')
#' 
#' 
#' # Size-based rarefaction and extrapolation sampling curves 
#' data(Brazil_rainforests)
#' data(Brazil_distM)
#' output4s = iNEXTbeta3D(data = Brazil_rainforests, diversity = 'FD', 
#'                        datatype = 'abundance', base = "size", nboot = 10, conf = 0.95, 
#'                        FDdistM = Brazil_distM, FDtype = 'AUC', FDcut_number = 30)
#' 
#' ggiNEXTbeta3D(output4s)
#' 
#' 
#' ## Taxonomic diversity for incidence data
#' # Coverage-based rarefaction and extrapolation sampling curves 
#' data(Second_growth_forests)
#' output5c = iNEXTbeta3D(data = Second_growth_forests, diversity = 'TD', datatype = 'incidence_raw', 
#'                        base = "coverage", level = NULL, nboot = 20, conf = 0.95)
#' 
#' ggiNEXTbeta3D(output5c, type = 'B')
#' ggiNEXTbeta3D(output5c, type = 'D')
#' 
#' 
#' # Size-based rarefaction and extrapolation sampling curves 
#' data(Second_growth_forests)
#' output5s = iNEXTbeta3D(data = Second_growth_forests, diversity = 'TD', datatype = 'incidence_raw', 
#'                        base = "size", level = NULL, nboot = 30, conf = 0.95)
#' 
#' ggiNEXTbeta3D(output5s)
#' 
#' 
#' 
#' @export
ggiNEXTbeta3D = function(output, type = 'B'){
  
  transp = 0.4
  cbPalette <- rev(c("#999999", "#E69F00", "#56B4E9", "#009E73", 
                     "#330066", "#CC79A7", "#0072B2", "#D55E00"))
  
  if ((length(output[[1]]) == 7 & type == "B") | (length(output[[1]]) == 2)) {
    
    if (unique(output[[1]]$gamma$Diversity) == 'TD') { ylab = "Taxonomic diversity" }
    if (unique(output[[1]]$gamma$Diversity) == 'PD') { ylab = "Phylogenetic diversity" }
    if (unique(output[[1]]$gamma$Diversity) == 'meanPD') { ylab = "Mean phylogenetic diversity" }
    if (unique(output[[1]]$gamma$Diversity) == 'FD_tau') { ylab = "Functional diversity (given tau)" }
    if (unique(output[[1]]$gamma$Diversity) == 'FD_AUC') { ylab = "Functional diversity (AUC)" }
  }
  
  if (length(output[[1]]) == 7 & type == "D") {
    
    if (unique(output[[1]]$gamma$Diversity) == 'TD') { ylab = "Taxonomic dissimilarity" }
    if (unique(output[[1]]$gamma$Diversity) == 'PD') { ylab = "Phylogenetic dissimilarity" }
    if (unique(output[[1]]$gamma$Diversity) == 'meanPD') { ylab = "Mean phylogenetic dissimilarity" }
    if (unique(output[[1]]$gamma$Diversity) == 'FD_tau') { ylab = "Functional dissimilarity (given tau)" }
    if (unique(output[[1]]$gamma$Diversity) == 'FD_AUC') { ylab = "Functional dissimilarity (AUC)" }
  }
  
  if (length(output[[1]]) == 7) {
    if (type == 'B'){
      
      gamma = lapply(output, function(y) {
        
        tmp = y[["gamma"]]
        
        tmp[tmp == 'Observed_SC(n, gamma)'] = tmp[tmp == 'Observed_SC(T, gamma)'] = 'Observed'
        tmp[tmp == 'Extrap_SC(2n, gamma)'] = tmp[tmp == 'Extrap_SC(2T, gamma)'] = 'Extrapolation'
        
        tmp[tmp == 'Observed_SC(n, alpha)'] = ifelse( unique(tmp[tmp$Method == 'Observed_SC(n, alpha)',"Size"] < unique(tmp[tmp$Method == "Observed", "Size"]) ),
                                                      "Rarefaction", "Extrapolation")
        tmp[tmp == 'Extrap_SC(2n, alpha)'] = ifelse( unique(tmp[tmp$Method == 'Extrap_SC(2n, alpha)',"Size"] < unique(tmp[tmp$Method == "Observed", "Size"]) ),
                                                     "Rarefaction", "Extrapolation")
        
        tmp[tmp == 'Observed_SC(T, alpha)'] = ifelse( unique(tmp[tmp$Method == 'Observed_SC(T, alpha)',"Size"] < unique(tmp[tmp$Method == "Observed", "Size"]) ),
                                                      "Rarefaction", "Extrapolation")
        tmp[tmp == 'Extrap_SC(2T, alpha)'] = ifelse( unique(tmp[tmp$Method == 'Extrap_SC(2T, alpha)',"Size"] < unique(tmp[tmp$Method == "Observed", "Size"]) ),
                                                     "Rarefaction", "Extrapolation")
        
        tmp
        
      }) %>% do.call(rbind,.) %>% rename("Estimate" = "Gamma") %>% mutate(div_type = "Gamma") %>% as_tibble()
      
      alpha = lapply(output, function(y) {
        
        tmp = y[["alpha"]]
        
        tmp[tmp == 'Observed_SC(n, alpha)'] = tmp[tmp == 'Observed_SC(T, alpha)'] = 'Observed'
        tmp[tmp == 'Extrap_SC(2n, alpha)'] = tmp[tmp == 'Extrap_SC(2T, alpha)'] = 'Extrapolation'
        
        tmp[tmp == 'Observed_SC(n, gamma)'] = ifelse( unique(tmp[tmp$Method == 'Observed_SC(n, gamma)',"Size"] < unique(tmp[tmp$Method == "Observed", "Size"]) ),
                                                      "Rarefaction", "Extrapolation")
        tmp[tmp == 'Extrap_SC(2n, gamma)'] = ifelse( unique(tmp[tmp$Method == 'Extrap_SC(2n, gamma)',"Size"] < unique(tmp[tmp$Method == "Observed", "Size"]) ),
                                                     "Rarefaction", "Extrapolation")
        
        tmp[tmp == 'Observed_SC(T, gamma)'] = ifelse( unique(tmp[tmp$Method == 'Observed_SC(T, gamma)',"Size"] < unique(tmp[tmp$Method == "Observed", "Size"]) ),
                                                      "Rarefaction", "Extrapolation")
        tmp[tmp == 'Extrap_SC(2T, gamma)'] = ifelse( unique(tmp[tmp$Method == 'Extrap_SC(2T, gamma)',"Size"] < unique(tmp[tmp$Method == "Observed", "Size"]) ),
                                                     "Rarefaction", "Extrapolation")
        
        tmp
        
      }) %>% do.call(rbind,.) %>% rename("Estimate" = "Alpha") %>% mutate(div_type = "Alpha") %>% as_tibble()
      
      beta =  lapply(output, function(y) {
        
        tmp = y[["beta"]]
        
        tmp[tmp == 'Observed_SC(n, alpha)'] = tmp[tmp == 'Observed_SC(T, alpha)'] = 'Observed'
        tmp[tmp == 'Extrap_SC(2n, alpha)'] = tmp[tmp == 'Extrap_SC(2T, alpha)'] = 'Extrapolation'
        
        tmp[tmp == 'Observed_SC(n, gamma)'] = ifelse( unique(tmp[tmp$Method == 'Observed_SC(n, gamma)',"Size"] < unique(tmp[tmp$Method == "Observed", "Size"]) ),
                                                      "Rarefaction", "Extrapolation")
        tmp[tmp == 'Extrap_SC(2n, gamma)'] = ifelse( unique(tmp[tmp$Method == 'Extrap_SC(2n, gamma)',"Size"] < unique(tmp[tmp$Method == "Observed", "Size"]) ),
                                                     "Rarefaction", "Extrapolation")
        
        tmp[tmp == 'Observed_SC(T, gamma)'] = ifelse( unique(tmp[tmp$Method == 'Observed_SC(T, gamma)',"Size"] < unique(tmp[tmp$Method == "Observed", "Size"]) ),
                                                      "Rarefaction", "Extrapolation")
        tmp[tmp == 'Extrap_SC(2T, gamma)'] = ifelse( unique(tmp[tmp$Method == 'Extrap_SC(2T, gamma)',"Size"] < unique(tmp[tmp$Method == "Observed", "Size"]) ),
                                                     "Rarefaction", "Extrapolation")
        
        tmp
        
      })  %>% do.call(rbind,.) %>% rename("Estimate" = "Beta") %>% mutate(div_type = "Beta")  %>% as_tibble()
      
      
      # # Dropping out the points extrapolated over double reference size
      # gamma1 = data.frame() ; alpha1 = data.frame() ; beta1 = data.frame()
      # 
      # for(i in 1:length(unique(gamma$Dataset))){
      #   
      #   Gamma <- gamma %>% filter(Dataset==unique(gamma$Dataset)[i]) ; ref_size = unique(Gamma[Gamma$Method=="Observed",]$Size)
      #   Gamma = Gamma %>% filter(!(Order.q==0 & round(Size)>2*ref_size))
      #   
      #   Alpha <- alpha %>% filter(Dataset==unique(gamma$Dataset)[i]) ; Alpha = Alpha %>% filter(!(Order.q==0 & round(Size)>2*ref_size))
      #   Beta <- beta %>% filter(Dataset==unique(gamma$Dataset)[i]) ; Beta = Beta %>% filter(!(Order.q==0 & round(Size)>2*ref_size))
      #   
      #   gamma1 = rbind(gamma1,Gamma) ; alpha1 = rbind(alpha1,Alpha) ; beta1 = rbind(beta1,Beta)
      #   
      # }
      # 
      # gamma = gamma1 ; alpha = alpha1 ; beta= beta1
      
      df = rbind(gamma, alpha, beta)
      for (i in unique(gamma$Order.q)) df$Order.q[df$Order.q == i] = paste0('q = ', i)
      df$div_type <- factor(df$div_type, levels = c("Gamma","Alpha","Beta"))
      
      id_obs = which(df$Method == 'Observed')
      
      if (length(id_obs) > 0) {
        
        for (i in 1:length(id_obs)) {
          
          new = df[id_obs[i],]
          new$SC = new$SC - 0.0001
          new$Method = 'Rarefaction'
          
          newe = df[id_obs[i],]
          newe$SC = newe$SC + 0.0001
          newe$Method = 'Extrapolation'
          
          df = rbind(df, new, newe)
          
        }
      }
      
      
    }
    
    if (type == 'D'){
      
      C = lapply(output, function(y) {
        
        tmp = y[["1-C"]]
        
        tmp[tmp == 'Observed_SC(n, alpha)'] = tmp[tmp == 'Observed_SC(T, alpha)'] = 'Observed'
        tmp[tmp == 'Extrap_SC(2n, alpha)'] = tmp[tmp == 'Extrap_SC(2T, alpha)'] = 'Extrapolation'
        
        tmp[tmp == 'Observed_SC(n, gamma)'] = ifelse( unique(tmp[tmp$Method == 'Observed_SC(n, gamma)',"Size"] < unique(tmp[tmp$Method == "Observed", "Size"]) ),
                                                      "Rarefaction", "Extrapolation")
        tmp[tmp == 'Extrap_SC(2n, gamma)'] = ifelse( unique(tmp[tmp$Method == 'Extrap_SC(2n, gamma)',"Size"] < unique(tmp[tmp$Method == "Observed", "Size"]) ),
                                                     "Rarefaction", "Extrapolation")
        
        tmp[tmp == 'Observed_SC(T, gamma)'] = ifelse( unique(tmp[tmp$Method == 'Observed_SC(T, gamma)',"Size"] < unique(tmp[tmp$Method == "Observed", "Size"]) ),
                                                      "Rarefaction", "Extrapolation")
        tmp[tmp == 'Extrap_SC(2T, gamma)'] = ifelse( unique(tmp[tmp$Method == 'Extrap_SC(2T, gamma)',"Size"] < unique(tmp[tmp$Method == "Observed", "Size"]) ),
                                                     "Rarefaction", "Extrapolation")
        
        tmp
        
      }) %>% do.call(rbind,.) %>% rename("Estimate" = "Dissimilarity") %>% mutate(div_type = "1-CqN") %>% as_tibble()
      
      U = lapply(output, function(y) {
        
        tmp = y[["1-U"]]
        
        tmp[tmp == 'Observed_SC(n, alpha)'] = tmp[tmp == 'Observed_SC(T, alpha)'] = 'Observed'
        tmp[tmp == 'Extrap_SC(2n, alpha)'] = tmp[tmp == 'Extrap_SC(2T, alpha)'] = 'Extrapolation'
        
        tmp[tmp == 'Observed_SC(n, gamma)'] = ifelse( unique(tmp[tmp$Method == 'Observed_SC(n, gamma)',"Size"] < unique(tmp[tmp$Method == "Observed", "Size"]) ),
                                                      "Rarefaction", "Extrapolation")
        tmp[tmp == 'Extrap_SC(2n, gamma)'] = ifelse( unique(tmp[tmp$Method == 'Extrap_SC(2n, gamma)',"Size"] < unique(tmp[tmp$Method == "Observed", "Size"]) ),
                                                     "Rarefaction", "Extrapolation")
        
        tmp[tmp == 'Observed_SC(T, gamma)'] = ifelse( unique(tmp[tmp$Method == 'Observed_SC(T, gamma)',"Size"] < unique(tmp[tmp$Method == "Observed", "Size"]) ),
                                                      "Rarefaction", "Extrapolation")
        tmp[tmp == 'Extrap_SC(2T, gamma)'] = ifelse( unique(tmp[tmp$Method == 'Extrap_SC(2T, gamma)',"Size"] < unique(tmp[tmp$Method == "Observed", "Size"]) ),
                                                     "Rarefaction", "Extrapolation")
        
        tmp
        
      }) %>% do.call(rbind,.) %>% rename("Estimate" = "Dissimilarity") %>% mutate(div_type = "1-UqN") %>% as_tibble()
      
      V = lapply(output, function(y) {
        
        tmp = y[["1-V"]]
        
        tmp[tmp == 'Observed_SC(n, alpha)'] = tmp[tmp == 'Observed_SC(T, alpha)'] = 'Observed'
        tmp[tmp == 'Extrap_SC(2n, alpha)'] = tmp[tmp == 'Extrap_SC(2T, alpha)'] = 'Extrapolation'
        
        tmp[tmp == 'Observed_SC(n, gamma)'] = ifelse( unique(tmp[tmp$Method == 'Observed_SC(n, gamma)',"Size"] < unique(tmp[tmp$Method == "Observed", "Size"]) ),
                                                      "Rarefaction", "Extrapolation")
        tmp[tmp == 'Extrap_SC(2n, gamma)'] = ifelse( unique(tmp[tmp$Method == 'Extrap_SC(2n, gamma)',"Size"] < unique(tmp[tmp$Method == "Observed", "Size"]) ),
                                                     "Rarefaction", "Extrapolation")
        
        tmp[tmp == 'Observed_SC(T, gamma)'] = ifelse( unique(tmp[tmp$Method == 'Observed_SC(T, gamma)',"Size"] < unique(tmp[tmp$Method == "Observed", "Size"]) ),
                                                      "Rarefaction", "Extrapolation")
        tmp[tmp == 'Extrap_SC(2T, gamma)'] = ifelse( unique(tmp[tmp$Method == 'Extrap_SC(2T, gamma)',"Size"] < unique(tmp[tmp$Method == "Observed", "Size"]) ),
                                                     "Rarefaction", "Extrapolation")
        
        tmp
        
      }) %>% do.call(rbind,.) %>% rename("Estimate" = "Dissimilarity") %>% mutate(div_type = "1-VqN") %>% as_tibble()
      
      S = lapply(output, function(y) {
        
        tmp = y[["1-S"]]
        
        tmp[tmp == 'Observed_SC(n, alpha)'] = tmp[tmp == 'Observed_SC(T, alpha)'] = 'Observed'
        tmp[tmp == 'Extrap_SC(2n, alpha)'] = tmp[tmp == 'Extrap_SC(2T, alpha)'] = 'Extrapolation'
        
        tmp[tmp == 'Observed_SC(n, gamma)'] = ifelse( unique(tmp[tmp$Method == 'Observed_SC(n, gamma)',"Size"] < unique(tmp[tmp$Method == "Observed", "Size"]) ),
                                                      "Rarefaction", "Extrapolation")
        tmp[tmp == 'Extrap_SC(2n, gamma)'] = ifelse( unique(tmp[tmp$Method == 'Extrap_SC(2n, gamma)',"Size"] < unique(tmp[tmp$Method == "Observed", "Size"]) ),
                                                     "Rarefaction", "Extrapolation")
        
        tmp[tmp == 'Observed_SC(T, gamma)'] = ifelse( unique(tmp[tmp$Method == 'Observed_SC(T, gamma)',"Size"] < unique(tmp[tmp$Method == "Observed", "Size"]) ),
                                                      "Rarefaction", "Extrapolation")
        tmp[tmp == 'Extrap_SC(2T, gamma)'] = ifelse( unique(tmp[tmp$Method == 'Extrap_SC(2T, gamma)',"Size"] < unique(tmp[tmp$Method == "Observed", "Size"]) ),
                                                     "Rarefaction", "Extrapolation")
        
        tmp
        
      }) %>% do.call(rbind,.) %>% rename("Estimate" = "Dissimilarity") %>% mutate(div_type = "1-SqN") %>% as_tibble()
      
      # C = C %>% filter(Method != 'Observed')
      # U = U %>% filter(Method != 'Observed')
      # V = V %>% filter(Method != 'Observed')
      # S = S %>% filter(Method != 'Observed')
     
      
      # # Dropping out the points extrapolated over double reference size
      # c1 = data.frame() ; u1 = data.frame() ; v1 = data.frame() ; s1 = data.frame()
      # 
      # for(i in 1:length(unique(C$Dataset))){
      #   
      #   CC <- C %>% filter(Dataset==unique(C$Dataset)[i]) ; ref_size = unique(CC[CC$Method=="Observed",]$Size)
      #   CC = CC %>% filter(!(Order.q==0 & round(Size)>2*ref_size))
      #   
      #   UU <- U %>% filter(Dataset==unique(C$Dataset)[i]) ; UU = UU %>% filter(!(Order.q==0 & round(Size)>2*ref_size))
      #   VV <- V %>% filter(Dataset==unique(C$Dataset)[i]) ; VV = VV %>% filter(!(Order.q==0 & round(Size)>2*ref_size))
      #   SS <- S %>% filter(Dataset==unique(C$Dataset)[i]) ; SS = SS %>% filter(!(Order.q==0 & round(Size)>2*ref_size))
      #   
      #   c1 = rbind(c1,CC) ; u1 = rbind(u1,UU) ; v1 = rbind(v1,VV) ; s1 = rbind(s1,SS)
      #   
      # }
      # 
      # C = c1 ; U = u1 ; V = v1 ; S = s1
      
      
      df = rbind(C, U, V, S)
      for (i in unique(C$Order.q)) df$Order.q[df$Order.q == i] = paste0('q = ', i)
      df$div_type <- factor(df$div_type, levels = c("1-CqN", "1-UqN", "1-VqN", "1-SqN"))
      
      id_obs = which(df$Method == 'Observed')
      
      if (length(id_obs) > 0) {
        
        for (i in 1:length(id_obs)) {
          
          new = df[id_obs[i],]
          new$SC = new$SC - 0.0001
          new$Method = 'Rarefaction'
          
          newe = df[id_obs[i],]
          newe$SC = newe$SC + 0.0001
          newe$Method = 'Extrapolation'
          
          df = rbind(df, new, newe)
          
        }
      }

      
    }
    
    lty = c(Rarefaction = "solid", Extrapolation = "dashed")
    df$Method = factor(df$Method, levels = c('Rarefaction', 'Extrapolation', 'Observed'))
    
    double_size = unique(df[df$Method == "Observed",]$Size)*2
    double_extrapolation = df %>% filter(Method == "Extrapolation" & round(Size) %in% double_size)
    
    fig = ggplot(data = df, aes(x = SC, y = Estimate, col = Dataset)) +
      geom_line(data = subset(df, Method != 'Observed'), aes(linetype = Method), size=1.1) + scale_linetype_manual(values = lty) +
      # geom_line(lty=2) + 
      geom_point(data = subset(df, Method == 'Observed' & div_type == "Gamma"), shape = 19, size = 3) + 
      geom_point(data = subset(df, Method == 'Observed' & div_type != "Gamma"), shape = 1, size = 3, stroke = 1.5)+
      geom_point(data = subset(double_extrapolation, div_type == "Gamma"), shape = 17, size = 3) + 
      geom_point(data = subset(double_extrapolation, div_type != "Gamma"), shape = 2, size = 3, stroke = 1.5) + 
      scale_colour_manual(values = cbPalette) + 
      scale_fill_manual(values = cbPalette) + 
      facet_grid(div_type ~ Order.q, scales = 'free') +
      theme_bw() + 
      theme(legend.position = "bottom", 
            legend.title = element_blank(),
            strip.text = element_text(size = 15, face = 'bold'),
            axis.title = element_text(hjust = 0.5, size = 15, face = 'bold'),
            axis.text.x = element_text(size = 12),
            axis.text.y = element_text(size = 12),
            legend.box = "vertical",
            legend.margin = margin(0, 0, 0, 0),
            legend.box.margin = margin(-10, -10, -5, -10),
            legend.text = element_text(size = 13),
            plot.margin = unit(c(5.5, 5.5, 5.5, 5.5), "pt")) +
      labs(x = 'Sample coverage', y = ylab) +
      guides(linetype = guide_legend(keywidth = 2.5))
    
    if (sum(is.na(df$LCL) == 0)) fig = fig + geom_ribbon(aes(ymin = LCL, ymax = UCL, fill = Dataset, col = NULL), alpha = transp)
    
  } else if (length(output[[1]]) == 2) {
    
    gamma = lapply(output, function(y) y[["gamma"]]) %>% do.call(rbind,.) %>% rename("Estimate" = "Gamma") %>% mutate(div_type = "Gamma") %>% as_tibble()
    alpha = lapply(output, function(y) y[["alpha"]]) %>% do.call(rbind,.) %>% rename("Estimate" = "Alpha") %>% mutate(div_type = "Alpha") %>% as_tibble()
    
    if ('nT' %in% colnames(gamma)) {
      
      xlab = 'Number of sampling units'
      colnames(gamma)[colnames(gamma) == 'nT'] = 'Size'
      colnames(alpha)[colnames(alpha) == 'nT'] = 'Size'
      
    } else xlab = 'Number of individuals'
    
    df = rbind(gamma, alpha)
    for (i in unique(gamma$Order.q)) df$Order.q[df$Order.q == i] = paste0('q = ', i)
    df$div_type <- factor(df$div_type, levels = c("Gamma","Alpha"))
    
    id_obs = which(df$Method == 'Observed')
    
    if (sum(id_obs) > 0) {
      for (i in 1:length(id_obs)) {
        
        new = df[id_obs[i],]
        new$Size = new$Size - 0.0001
        new$Method = 'Rarefaction'
        
        newe = df[id_obs[i],]
        newe$Size = newe$Size + 0.0001
        newe$Method = 'Extrapolation'
        
        df = rbind(df, new, newe)
        
      }
    }
    
    lty = c(Rarefaction = "solid", Extrapolation = "dashed")
    df$Method = factor(df$Method, levels = c('Rarefaction', 'Extrapolation', 'Observed'))
    
    double_size = unique(df[df$Method == "Observed",]$Size)*2
    double_extrapolation = df %>% filter(Method == "Extrapolation" & round(Size) %in% double_size)
    
    fig = ggplot(data = df, aes(x = Size, y = Estimate, col = Dataset)) +
      geom_line(data = subset(df, Method != 'Observed'), aes(linetype = Method), size=1.1) + scale_linetype_manual(values = lty) +
      geom_point(data = subset(df, Method == 'Observed' & div_type == "Gamma"), shape = 19, size = 3) + 
      geom_point(data = subset(df, Method == 'Observed' & div_type != "Gamma"), shape = 1, size = 3, stroke = 1.5) +
      geom_point(data = subset(double_extrapolation, div_type == "Gamma"), shape = 17, size = 3) + 
      geom_point(data = subset(double_extrapolation, div_type != "Gamma"), shape = 2, size = 3, stroke = 1.5) + 
      scale_colour_manual(values = cbPalette) + 
      scale_fill_manual(values = cbPalette) + 
      facet_grid(div_type ~ Order.q, scales = 'free') +
      theme_bw() + 
      theme(legend.position = "bottom", 
            legend.title = element_blank(),
            strip.text = element_text(size = 15, face = 'bold'),
            axis.title = element_text(hjust = 0.5, size = 15, face = 'bold'),
            axis.text.x = element_text(size = 12),
            axis.text.y = element_text(size = 12),
            legend.box = "vertical",
            legend.margin = margin(0, 0, 0, 0),
            legend.box.margin = margin(-10, -10, -5, -10),
            legend.text = element_text(size = 13),
            plot.margin = unit(c(5.5, 5.5, 5.5, 5.5), "pt")) +
      labs(x = xlab, y = ylab) +
      guides(linetype = guide_legend(keywidth = 2.5))
    
    if (sum(is.na(df$LCL) == 0)) fig = fig + geom_ribbon(aes(ymin = LCL, ymax = UCL, fill = Dataset, col = NULL), alpha = transp)
  }
  
  return(fig)
}


coverage_to_size = function(x, C, datatype = 'abundance'){
  
  if (datatype == 'abundance'){
    
    n <- sum(x)
    refC <- iNEXT.3D:::Coverage(x, 'abundance', n)
    f <- function(m, C) abs(iNEXT.3D:::Coverage(x, 'abundance', m) - C)
    if (refC == C) {
      mm = n
    } else if (refC > C) {
      opt <- optimize(f, C = C, lower = 0, upper = sum(x))
      mm <- opt$minimum
      # mm <- round(mm)
    } else if (refC < C) {
      f1 <- sum(x == 1)
      f2 <- sum(x == 2)
      if (f1 > 0 & f2 > 0) {
        A <- (n - 1) * f1/((n - 1) * f1 + 2 * f2)
      }
      if (f1 > 1 & f2 == 0) {
        A <- (n - 1) * (f1 - 1)/((n - 1) * (f1 - 1) + 2)
      }
      if (f1 == 1 & f2 == 0) {
        A <- 1
      }
      if (f1 == 0 & f2 == 0) {
        A <- 1
      }
      if (f1 == 0 & f2 > 0) {
        A <- 1
      }
      mm <- (log(n/f1) + log(1 - C))/log(A) - 1
      if (is.nan(mm) == TRUE) mm = Inf
      mm <- n + mm
      # mm <- round(mm)
    }
    
  } else {
    
    m <- NULL
    n <- max(x)
    refC <- iNEXT.3D:::Coverage(x, 'incidence_freq', n)
    f <- function(m, C) abs(iNEXT.3D:::Coverage(x, 'incidence_freq', m) - C)
    if (refC == C) {
      mm = n
    } else if (refC > C) {
      opt <- optimize(f, C = C, lower = 0, upper = max(x))
      mm <- opt$minimum
      # mm <- round(mm)
    } else if (refC < C) {
      f1 <- sum(x == 1)
      f2 <- sum(x == 2)
      U <- sum(x) - max(x)
      if (f1 > 0 & f2 > 0) {
        A <- (n - 1) * f1/((n - 1) * f1 + 2 * f2)
      }
      if (f1 > 1 & f2 == 0) {
        A <- (n - 1) * (f1 - 1)/((n - 1) * (f1 - 1) + 2)
      }
      if (f1 == 1 & f2 == 0) {
        A <- 1
      }
      if (f1 == 0 & f2 == 0) {
        A <- 1
      }
      if (f1 == 0 & f2 > 0) {
        A <- 1
      }
      mm <- (log(U/f1) + log(1 - C))/log(A) - 1
      if (is.nan(mm) == TRUE) mm = Inf
      mm <- n + mm
      # mm <- round(mm)
    }
    
  }
  
  return(mm)
}

bootstrap_population_multiple_assemblage = function(data, data_gamma, datatype){
  
  if (datatype == 'abundance'){
    
    S_obs = sum(data_gamma > 0)
    n = sum(data_gamma)
    f1 = sum(data_gamma == 1)
    f2 = sum(data_gamma == 2)
    f0_hat = ifelse(f2 == 0, (n - 1)/n * f1 * (f1 - 1)/2, (n - 1)/n * f1^2/2/f2) %>% ceiling()
    
    output = apply(data, 2, function(x){
      
      p_i_hat = iNEXT.3D:::EstiBootComm.Ind(Spec = x)
      
      if(length(p_i_hat) != length(x)){
        
        p_i_hat_unobs = p_i_hat[(length(x)+1):length(p_i_hat)]
        p_i_hat_obs = p_i_hat[1:length(x)]
        p_i_hat = c(p_i_hat_obs, rep(0, f0_hat))
        candidate = which(p_i_hat==0)
        
        chosen = sample(x = candidate, size = min(length(p_i_hat_unobs), length(candidate)), replace = F)
        p_i_hat[chosen] = (1-sum(p_i_hat))/length(chosen)
        
        p_i_hat
        
      } else {
        
        p_i_hat = c(p_i_hat, rep(0, f0_hat))
        p_i_hat
        
      }
    })
    
  }
  
  if (datatype == 'incidence'){
    
    S_obs = sum(data_gamma > 0)
    t = data_gamma[1]
    Q1 = sum(data_gamma == 1)
    Q2 = sum(data_gamma == 2)
    Q0_hat = if ( Q2 == 0 ){( (t-1)/t ) * ( Q1*(Q1-1)/2 )} else {( (t-1)/t ) * ( (Q1^2) / (2*Q2) )} %>% ceiling
    
    output = apply(data, 2, function(x){
      
      pi_i_hat = iNEXT.3D:::EstiBootComm.Sam(Spec = x)
      
      if(length(pi_i_hat) != (length(x) - 1)){
        
        pi_i_hat_unobs = pi_i_hat[length(x):length(pi_i_hat)]
        pi_i_hat_obs = pi_i_hat[1:(length(x)-1)]
        pi_i_hat = c(pi_i_hat_obs, rep(0, Q0_hat))
        candidate = which(pi_i_hat == 0)
        
        chosen = sample(x = candidate, size = min(length(pi_i_hat_unobs), length(candidate)), replace = F)
        pi_i_hat[chosen] = unique(pi_i_hat_unobs)
        
        pi_i_hat
        
      } else {
        
        pi_i_hat = c(pi_i_hat, rep(0, Q0_hat))
        pi_i_hat
        
      }
    })
    
  }
  
  return(output)
  
}

Bootstrap_distance_matrix = function(data, distance_matrix, f0.hat, datatype){
  
  if (datatype == "incidence_freq") {
    n = data[1]
    X = data[-1]
    u = sum(data)
  } else if (datatype == "abundance") {
    n = sum(data)
    X = data
  }
  
  # n = sum(data)
  distance = as.matrix(distance_matrix)
  dij = distance
  # X = data
  
  F.1 <- sum(dij[, X==1]) ; F.2 <- sum(dij[, X==2])
  F11 <- sum(dij[X==1, X==1]) ; F22 <- sum(dij[X==2, X==2])
  
  if (datatype == "abundance") {
    F.0hat <- ifelse(F.2 > 0, ((n-1)/n) * (F.1^2/(2 * F.2)), ((n-1)/n)*(F.1*(F.1-0.01)/(2)))
    F00hat <- ifelse(F22 > 0, ((n-2)* (n-3)* (F11^2)/(4* n* (n-1)* F22)), ((n-2)* (n-3)* (F11*(F11-0.01))/(4 *n * (n-1))) )
  } else if (datatype == "incidence_freq") {
    F.0hat <- ifelse(F.2 > 0, ((n-1)/n) * (F.1^2/(2 * F.2)), ((n-1)/n)*(F.1*(F.1-0.01)/(2)))
    F00hat <- ifelse(F22 > 0, ((n-1)^2 * (F11^2)/(4* n* n* F22)), ((n-1)* (n-1)* (F11*(F11-0.01))/(4 *n * n)) )
  }
  
  if (f0.hat == 0) {
    d = dij
  } else if (f0.hat == 1) {
    d.0bar <- matrix(rep(F.0hat/length(X)/f0.hat, length(X)*f0.hat), length(X), f0.hat)
    
    d00 = matrix(0, f0.hat, f0.hat)
    d <- cbind(dij, d.0bar )
    aa <- cbind(t(d.0bar), d00 )
    d <- rbind(d, aa)
    diag(d) = 0
  } else {
    d.0bar <- matrix(rep(F.0hat/length(X)/f0.hat, length(X)*f0.hat), length(X), f0.hat)
    
    fo.num = (f0.hat * (f0.hat-1) )/2
    d00 = matrix(0, f0.hat, f0.hat)
    d00[upper.tri(d00)] = (F00hat/2)/fo.num
    d00 <- pmax(d00, t(d00))###signmatrix
    d <- cbind(dij, d.0bar )
    aa <- cbind(t(d.0bar), d00 )
    d <- rbind(d, aa)
    diag(d) = 0
  }
  
  return(d)
  
}

FD.m.est_0 = function (ai_vi, m, q, nT) {
  EFD = function(m, qs, obs, asy, beta, av) {
    m = m - nT
    out <- sapply(1:length(qs), function(i) {
      if (qs[i] != 2) {
        obs[i] + (asy[i] - obs[i]) * (1 - (1 - beta[i])^m)
      }
      else if (qs[i] == 2) {
        V_bar^2/sum((av[, 2]) * ((1/(nT + m)) * (av[, 1]/nT) + ((nT + m - 1)/(nT + m)) * (av[, 1] * (av[, 1] - 1)/(nT * (nT - 1)))))
      }
    })
    return(out)
  }
  V_bar <- sum(ai_vi$ai[, 1] * ai_vi$vi[, 1])/nT
  asy <- iNEXT.3D:::FD_est(ai_vi, q, nT)$est
  obs <- iNEXT.3D:::FD_mle(ai_vi, q)
  out <- sapply(1:ncol(ai_vi$ai), function(i) {
    ai <- ai_vi$ai[, i]
    ai[ai < 1] <- 1
    av = cbind(ai = round(ai), vi = ai_vi$vi[, i])
    RFD_m = iNEXT.3D:::RFD(av, nT, nT - 1, q, V_bar)
    beta <- rep(0, length(q))
    asy_i <- asy[, i]
    obs_i <- obs[, i]
    asy_i <- sapply(1:length(q), function(j) {
      max(asy_i[j], obs_i[j], RFD_m[j])
    })
    
    obs_i <- sapply(1:length(q), function(j) {
      max(RFD_m[j], obs_i[j])
    })
    
    beta0plus <- which(asy_i != obs_i)
    beta[beta0plus] <- (obs_i[beta0plus] - RFD_m[beta0plus])/(asy_i[beta0plus] - RFD_m[beta0plus])
    
    if (sum(m < nT) != 0) {
      int.m = sort(unique(c(floor(m[m < nT]), ceiling(m[m < nT]))))
      mRFD = rbind(int.m, sapply(int.m, function(k) iNEXT.3D:::RFD(av, nT, k, q, V_bar)))
      
      if (0 %in% int.m) mRFD[,mRFD[1,] == 0] = 0
    }
    
    sapply(m, function(mm) {
      if (mm < nT) {
        if (mm == round(mm)) {
          mRFD[-1, mRFD[1, ] == mm]
        } else {
          (ceiling(mm) - mm) * mRFD[-1, mRFD[1, ] == floor(mm)] + (mm - floor(mm)) * mRFD[-1, mRFD[1, ] == ceiling(mm)]
        }
      }
      else if (mm == nT) {
        obs_i
      }
      else if (mm == Inf) {
        asy_i
      }
      else {
        EFD(m = mm, qs = q, obs = obs_i, asy = asy_i, beta = beta, av = av)
      }
    }) %>% t %>% as.numeric
  })
  matrix(out, ncol = ncol(ai_vi$ai))
}



#' Data information for reference samples 
#' 
#' \code{DataInfobeta3D}: Provides basic data information for (1) the reference sample in each assemblage, 
#' (2) the gamma reference sample in the pooled assemblage, and (3) the alpha reference sample in the
#'  joint assemblage for TD, mean-PD and FD. 
#' 
#' @param data (a) For \code{datatype = "abundance"}, species abundance data for a single dataset can be input as a \code{matrix/data.frame} (species-by-assemblage); data for multiple datasets can be input as a \code{list} of \code{matrices/data.frames}, with each matrix representing a species-by-assemblage abundance matrix for one of the datasets.\cr
#' (b) For \code{datatype = "incidence_raw"}, data for a single dataset with N assemblages can be input as a \code{list} of \code{matrices/data.frames}, with each matrix representing a species-by-sampling-unit incidence matrix for one of the assemblages; data for multiple datasets can be input as multiple lists.
#' @param diversity selection of diversity type: \code{'TD'} = Taxonomic diversity, \code{'PD'} = Phylogenetic diversity, and \code{'FD'} = Functional diversity.
#' @param datatype data type of input data: individual-based abundance data (\code{datatype = "abundance"}) or species by sampling-units incidence matrix (\code{datatype = "incidence_raw"}) with all entries being \code{0} (non-detection) or \code{1} (detection).
#' @param PDtree (required argument only for \code{diversity = "PD"}), a phylogenetic tree in Newick format for all observed species in the pooled dataset. 
#' @param PDreftime (argument only for \code{diversity = "PD"}), a numerical value specifying reference time for PD. Default is \code{PDreftime = NULL} (i.e., the age of the root of PDtree).  
#' @param FDdistM (required argument only for \code{diversity = "FD"}), a species pairwise distance matrix for all species in the pooled dataset. 
#' @param FDtype (argument only for \code{diversity = "FD"}), select FD type: \code{FDtype = "tau_value"} for FD under a specified threshold value, or \code{FDtype = "AUC"} (area under the curve of tau-profile) for an overall FD which integrates all threshold values between zero and one. Default is \code{FDtype = "AUC"}.  
#' @param FDtau (argument only for \code{diversity = "FD"} and \code{FDtype = "tau_value"}), a numerical value between 0 and 1 specifying the tau value (threshold level). If \code{FDtype = NULL} (default), then threshold is set to be the mean distance between any two individuals randomly selected from the pooled dataset (i.e., quadratic entropy). 
#' 
#' @return a data.frame including basic data information.\cr\cr 
#' For abundance data, basic information shared by TD, mean-PD and FD
#'  includes dataset name (Dataset), individual/pooled/joint assemblage (Assemblage),
#' sample size (n), observed species richness (S.obs), sample coverage estimates of the reference sample (SC(n)), 
#' sample coverage estimate for twice the reference sample size (SC(2n)). Other additional information is given below.\cr\cr
#' (1) TD: the first five species abundance (f1-f5).\cr\cr
#' (2) Mean-PD: the the observed total branch length in the phylogenetic tree (PD.obs), 
#' the number of singletons (f1*) and doubletons (f2*) in the node/branch abundance set, as well as the total branch length 
#' of those singletons (g1) and of those doubletons (g2), and the reference time (Reftime).\cr\cr
#' (3) FD (\code{FDtype = "AUC"}): the minimum distance among all non-diagonal elements in the distance matrix (dmin), the mean distance
#' (dmean), and the maximum distance (dmax) in the distance matrix.\cr \cr
#' (4) FD (\code{FDtype = "tau_value"}): the number of singletons (a1*) and of doubletons (a2*) among the functionally indistinct
#'  set at the specified threshold level 'Tau', as well as the total contribution of singletons (h1) and of doubletons (h2)
#'   at the specified threshold level 'Tau'.\cr\cr
#'  
#'  For incidence data, the basic information for TD includes dataset name (Dataset), individual/pooled/joint assemblage 
#'  (Assemblage), number of sampling units (T), total number of incidences (U), observed species richness (S.obs), 
#'  sample coverage estimates of the reference sample (SC(T)), sample coverage estimate for twice the reference sample size
#'  (SC(2T)), as well as the first species incidence frequency counts (Q1-Q5). For mean-PD and FD, output is similar to that
#'  for abundance data.   
#' 
#' @examples
#' ## Taxonomic diversity for abundance data
#' data(Brazil_rainforests)
#' output1 = DataInfobeta3D(data = Brazil_rainforests, diversity = 'TD', datatype = 'abundance')
#' output1
#' 
#' 
#' ## Mean phylogenetic diversity for abundance data
#' data(Brazil_rainforests)
#' data(Brazil_tree)
#' output2 = DataInfobeta3D(data = Brazil_rainforests, diversity = 'PD', 
#'                          datatype = 'abundance', PDtree = Brazil_tree, PDreftime = NULL)
#' output2
#' 
#' 
#' ## Functional diversity for abundance data under a specified threshold level
#' data(Brazil_rainforests)
#' data(Brazil_distM)
#' output3 = DataInfobeta3D(data = Brazil_rainforests, diversity = 'FD', 
#'                          datatype = 'abundance', FDdistM = Brazil_distM, FDtype = 'tau_value', FDtau = NULL)
#' output3
#' 
#' 
#' ## Functional diversity for abundance data when all threshold levels from 0 to 1 are considered
#' data(Brazil_rainforests)
#' data(Brazil_distM)
#' output4 = DataInfobeta3D(data = Brazil_rainforests, diversity = 'FD', 
#'                          datatype = 'abundance', FDdistM = Brazil_distM, FDtype = 'AUC')
#' output4
#' 
#' 
#' ## Taxonomic diversity for incidence data
#' data(Second_growth_forests)
#' output5 = DataInfobeta3D(data = Second_growth_forests, diversity = 'TD', datatype = 'incidence_raw')
#' output5
#' 
#' 
#' @export
DataInfobeta3D = function(data, diversity = 'TD', datatype = 'abundance', 
                          PDtree = NULL, PDreftime = NULL, FDdistM = NULL, FDtype = 'AUC', FDtau = NULL) {
  
  ## Check parameter setting
  if (is.na(pmatch(diversity, c("TD", "PD", "FD")))) stop("invalid diversity")
  
  if (datatype == "incidence" | datatype == "incidence_freq") stop('Please try datatype = "incidence_raw".')  
  if(is.na(pmatch(datatype, c("abundance", "incidence_raw"))))
    stop("invalid datatype")
  
  if (! (is.null(PDreftime) | inherits(PDreftime, "numeric")))
    stop("invalid class of reference time, PDreftime should be a postive value of numeric object.", call. = FALSE)
  if (length(PDreftime) > 1)
    stop("PDreftime can only accept a value instead of a vector.", call. = FALSE)
  
  if (FDtype == "tau_values") stop('Please try FDtype = "tau_value".')  
  if(is.na(pmatch(FDtype, c("AUC", "tau_value"))))
    stop("Incorrect type of functional diversity type, please use either 'AUC' or 'tau_value'", call. = FALSE)
  
  if (! (is.null(FDtau) | inherits(FDtau, "numeric")))
    stop("invalid class of tau value, FDtau should be a postive value between zero and one.", call. = FALSE)
  if (length(FDtau) > 1)
    stop("FDtau only accept a value instead of a vector.", call. = FALSE)
  ##
  
  if (is.null(names(data))) names(data) = paste0("Dataset_", 1:length(data))
  
  if (datatype == "abundance") data = lapply(data, as.matrix) else if (datatype == "incidence_raw")
    data = lapply(data, function(x) lapply(x, as.matrix))
  
  if (diversity == "TD") {
    
    if (datatype == "abundance") {
      
      Dat = lapply(data, function(x) list("Pooled assemblage" = rowSums(x),
                                          "Joint assemblage" = as.vector(x)))
      
      output = lapply(1:length(Dat), function(i) {
        
        if (is.null(colnames(data[[i]]))) colnames(data[[i]]) = paste('Assemblage_', 1:ncol(data[[i]]), sep = "")
        
        multiple = DataInfo3D(Dat[[i]], datatype = "abundance")
        single = DataInfo3D(data[[i]], datatype = "abundance")
        
        
        return(rbind(single, multiple) %>% cbind(Dataset = names(data)[i],.) %>% .[,1:11])
        
        }) %>% do.call(rbind,.)
    }
    
    if (datatype == "incidence_raw") {
      
      Dat = lapply(data, function(x) {
        
        data_gamma = Reduce('+', x)
        data_gamma[data_gamma > 1] = 1
        
        data_alpha = do.call(rbind, x)
        
        list("Pooled assemblage" = c(ncol(data_gamma), as.vector(rowSums(data_gamma))),
             "Joint assemblage" = c(ncol(data_alpha), as.vector(rowSums(data_alpha))))
      })
      
      output = lapply(1:length(Dat), function(i) {
        
        if (is.null(names(data[[i]]))) names(data[[i]]) = paste('Assemblage_', 1:length(data[[i]]), sep = "")
        
        multiple = DataInfo3D(Dat[[i]], datatype = "incidence_freq")
        single = DataInfo3D(data[[i]], datatype = "incidence_raw")
        
        return(rbind(single, multiple) %>% cbind(Dataset = names(data)[i],.) %>% .[,1:12])
        
        }) %>% do.call(rbind,.)
    }
    
  }
  
  if (diversity == "PD") {
    
    if(is.null(PDreftime)) PDreftime = get.rooted.tree.height(PDtree) else if (PDreftime <= 0) { 
      stop("Reference time must be greater than 0. Use NULL to set it to pooled tree height.", call. = FALSE)
      }
    
    
    if (datatype == "abundance") {
      
      output = lapply(1:length(data), function(i) {
        
        x = data[[i]]
        if (is.null(colnames(x))) colnames(x) = paste('Assemblage_', 1:ncol(x), sep = "")
        
        data_gamma = rowSums(x)
        data_gamma = data_gamma[data_gamma > 0]
        
        aL_table_gamma = iNEXT.3D:::phyBranchAL_Abu(phylo = PDtree, data = data_gamma, rootExtend = T, refT = PDreftime)
        aL_table_gamma$treeNabu$branch.length = aL_table_gamma$BLbyT[,1]
        aL_table_gamma = aL_table_gamma$treeNabu %>% select(branch.abun, branch.length, tgroup)
        
        N = ncol(x)
        aL_table_alpha = c()
        
        for (j in 1:N) {
          
          y = x[x[,j] > 0, j]
          names(y) = rownames(x)[x[,j]>0]
          
          aL_table = iNEXT.3D:::phyBranchAL_Abu(phylo = PDtree, data = y, rootExtend = T, refT = PDreftime)
          aL_table$treeNabu$branch.length = aL_table$BLbyT[,1]
          aL_table = aL_table$treeNabu %>% select(branch.abun, branch.length, tgroup)
          
          aL_table_alpha = rbind(aL_table_alpha, aL_table)
          
        }
        
        output <- sapply(list(aL_table_gamma, aL_table_alpha), function(y){
          
          ai = y$branch.abun
          Li = y$branch.length
          I1 <- which(ai == 1 & Li > 0)
          I2 <- which(ai == 2 & Li > 0)
          S.obs = y[y$branch.abun > 0,] %>% filter(tgroup == "Tip") %>% nrow
          PD_obs <- sum(Li)
          f1 <- length(I1)
          f2 <- length(I2)
          g1 <- sum(Li[I1])
          g2 <- sum(Li[I2])
          c(S.obs, PD_obs, f1, f2, g1, g2)
          
        }) %>% t()
        
        Chat = sapply(list(rowSums(x), as.vector(x)), function(y) iNEXT.3D:::Coverage(y, "abundance", sum(y)))
        mul_C2n = sapply(list(rowSums(x), as.vector(x)), function(y) iNEXT.3D:::Coverage(y, "abundance", 2*sum(y)))
        
        multiple = tibble('Assemblage' = c("Pooled assemblage", "Joint assemblage"), 
                          'n' = sum(x), 'S.obs' = output[,1], 'SC(n)' = Chat, 'SC(2n)' = mul_C2n, 'PD.obs' = output[,2],
                          'f1*' = output[,3], 'f2*' = output[,4], 'g1' = output[,5], 'g2' = output[,6],
                          'Reftime' = PDreftime)
        
        single = DataInfo3D(x, diversity = "PD", datatype = "abundance", PDtree = PDtree, PDreftime = PDreftime)
        
        return(rbind(single, multiple) %>% cbind(Dataset = names(data)[i],.))
        
      }) %>% do.call(rbind,.)
      }
    
    if (datatype == "incidence_raw") {
      
      output = lapply(1:length(data), function(i) {
        
        x = data[[i]]
        if (is.null(names(x))) names(x) = paste('Assemblage_', 1:length(x), sep = "")
        
        data_gamma = Reduce('+', x)
        data_gamma[data_gamma > 1] = 1
        
        aL_table_gamma = iNEXT.3D:::phyBranchAL_Inc(phylo = PDtree, data = as.matrix(data_gamma), datatype = "incidence_raw", refT = PDreftime, rootExtend = T)
        aL_table_gamma$treeNabu$branch.length = aL_table_gamma$BLbyT[,1]
        aL_table_gamma = aL_table_gamma$treeNabu %>% select(branch.abun, branch.length, tgroup)
        
        N = length(x)
        aL_table_alpha = c()
        
        for (j in 1:N) {
          
          aL_table = iNEXT.3D:::phyBranchAL_Inc(phylo = PDtree, data = x[[j]], datatype = "incidence_raw", refT = PDreftime, rootExtend = T)
          aL_table$treeNabu$branch.length = aL_table$BLbyT[,1]
          aL_table = aL_table$treeNabu %>% select(branch.abun, branch.length, tgroup)
          
          aL_table_alpha = rbind(aL_table_alpha, aL_table)
          
        }
        
        output <- sapply(list(aL_table_gamma, aL_table_alpha), function(y){
          
          ai = y$branch.abun
          Li = y$branch.length
          I1 <- which(ai == 1 & Li > 0)
          I2 <- which(ai == 2 & Li > 0)
          S.obs = y[y$branch.abun > 0,] %>% filter(tgroup == "Tip") %>% nrow
          PD_obs <- sum(Li)
          f1 <- length(I1)
          f2 <- length(I2)
          g1 <- sum(Li[I1])
          g2 <- sum(Li[I2])
          c(S.obs, PD_obs, f1, f2, g1, g2)
          
          }) %>% t()
        
        Chat = sapply(list(data_gamma, do.call(rbind, x)), function(y) iNEXT.3D:::Coverage(y, "incidence_raw", ncol(y)))
        mul_C2T = sapply(list(data_gamma, do.call(rbind, x)), function(y) iNEXT.3D:::Coverage(y, "incidence_raw", 2*ncol(y)))
        multiple = tibble('Assemblage' = c("Pooled assemblage", "Joint assemblage"), 
                          'T' = ncol(x[[1]]), 'U' = c(sum(data_gamma), sum(sapply(x, rowSums))), 
                          'S.obs' = output[,1], 'SC(T)' = Chat, 'SC(2T)' = mul_C2T, 'PD.obs' = output[,2],
                          'Q1*' = output[,3], 'Q2*' = output[,4], 'R1' = output[,5], 'R2' = output[,6],
                          'Reftime' = PDreftime)
        
        single = DataInfo3D(x, diversity = "PD", datatype = "incidence_raw", PDtree = PDtree, PDreftime = PDreftime)
        
        return(rbind(single, multiple) %>% cbind(Dataset = names(data)[i],.))
        
      }) %>% do.call(rbind,.)
      
    }
    
    rownames(output) = NULL
  }
  
  if (diversity == "FD" & FDtype == "tau_value") {
    
    FDdistM = as.matrix(FDdistM)
    
    if (is.null(FDtau) == T) {
      
      if (datatype == 'abundance') {
        
        tmp <- lapply(data, rowSums)
        tmp = lapply(tmp, function(i) data.frame('value' = i) %>% rownames_to_column(var = "Species"))
        pdata = tmp[[1]]
        for(i in 2:length(tmp)){
          pdata = full_join(pdata, tmp[[i]], by = "Species")
        }
        pdata[is.na(pdata)] = 0
        pdata = pdata %>% column_to_rownames("Species")
        pdata = rowSums(pdata)
        
        order_sp <- match(names(pdata),rownames(FDdistM))
        FDdistM <- FDdistM[order_sp,order_sp]
        pdata <- matrix(pdata/sum(pdata), ncol = 1)
        
      } else if (datatype == 'incidence_raw') {
        
        tmp <- lapply(data, function(x) {tmp = Reduce('+', x); tmp[tmp > 1] = 1; rowSums(tmp) })
        tmp = lapply(tmp, function(i) data.frame('value' = i) %>% rownames_to_column(var = "Species"))
        pdata = tmp[[1]]
        for(i in 2:length(tmp)){
          pdata = full_join(pdata, tmp[[i]], by = "Species")
        }
        pdata[is.na(pdata)] = 0
        pdata = pdata %>% column_to_rownames("Species")
        pdata = rowSums(pdata)
        
        order_sp <- match(names(pdata),rownames(FDdistM))
        FDdistM <- FDdistM[order_sp,order_sp]
        pdata <- matrix(pdata/sum(pdata), ncol = 1)
        
      }
      
      FDtau <- sum ( (pdata %*% t(pdata) ) * FDdistM) # dmean
    }
    
    
    if (datatype == "abundance") {
      
      output = lapply(1:length(data), function(i) {
        
        x = data[[i]]
        if (is.null(colnames(x))) colnames(x) = paste('Assemblage_', 1:ncol(x), sep = "")
        
        dij = FDdistM
        dij = dij[rownames(dij) %in% rownames(x), colnames(dij) %in% rownames(x)]
        order_sp <- match(rownames(x), rownames(dij))
        dij <- dij[order_sp, order_sp]
        
        dij = dij[rowSums(x)>0, rowSums(x)>0]
        x = x[rowSums(x)>0,]
        
        
        if (FDtau == 0) {
          dij[dij > 0] <- 1
          aik = (1 - dij/1) %*% as.matrix(x)
        } else {
          dij[which(dij > FDtau, arr.ind = T)] = FDtau
          aik = (1 - dij/FDtau) %*% as.matrix(x)
        }
        
        positive_id = rowSums(aik) > 0
        
        
        gamma_x = rowSums(x)[positive_id]
        gamma_a = rowSums(aik)[positive_id]
        gamma_v = gamma_x/gamma_a
        
        alpha_a = as.vector(aik)
        alpha_v = rep(gamma_v, ncol(x))
        
        Chat = sapply(list(rowSums(x), as.vector(x)), function(y) iNEXT.3D:::Coverage(y, "abundance", sum(y)))
        mul_C2n = sapply(list(rowSums(x), as.vector(x)), function(y) iNEXT.3D:::Coverage(y, "abundance", 2*sum(y)))
        
        multiple = tibble('Assemblage' = c("Pooled assemblage", "Joint assemblage"), 
                          'n' = sum(x), 
                          'S.obs' = c(sum(rowSums(x) > 0), sum(as.vector(x) > 0)), 
                          'SC(n)' = Chat, 'SC(2n)' = mul_C2n,
                          'a1*' = c(sum(round(gamma_a) == 1), sum(round(alpha_a) == 1)), 
                          'a2*' = c(sum(round(gamma_a) == 2), sum(round(alpha_a) == 2)), 
                          'h1' = c(sum(gamma_v[round(gamma_a) == 1]), sum(alpha_v[round(alpha_a) == 1])), 
                          'h2' = c(sum(gamma_v[round(gamma_a) == 2]), sum(alpha_v[round(alpha_a) == 2])),
                          'Tau' = FDtau)
        
        single = DataInfo3D(x, diversity = "FD", datatype = "abundance", FDtype = "tau_values", FDtau = FDtau, FDdistM = FDdistM)
        
        return(rbind(single, multiple) %>% cbind(Dataset = names(data)[i],.))
        
      }) %>% do.call(rbind,.)
    }
    
    if (datatype == "incidence_raw") {
      
      output = lapply(1:length(data), function(i) {
        
        x = data[[i]]
        if (is.null(names(x))) names(x) = paste('Assemblage_', 1:length(x), sep = "")
        
        nT = ncol(x[[1]])
        
        data_gamma = Reduce('+', x)
        data_gamma[data_gamma > 1] = 1
        data_gamma_freq = rowSums(data_gamma)
        
        data_2D = sapply(x, rowSums)
        
        ##
        dij = FDdistM
        dij = dij[rownames(dij) %in% names(data_gamma_freq), colnames(dij) %in% names(data_gamma_freq)]
        order_sp <- match(names(data_gamma_freq), rownames(dij))
        dij <- dij[order_sp, order_sp]
        
        dij = dij[data_gamma_freq > 0, data_gamma_freq > 0]
        data_gamma_freq = data_gamma_freq[data_gamma_freq > 0]
        
        
        if (FDtau == 0) {
          dij[dij > 0] <- 1
          gamma_a = (1 - dij/1) %*% as.matrix(data_gamma_freq)
        } else {
          dij[which(dij > FDtau, arr.ind = T)] = FDtau
          gamma_a = (1 - dij/FDtau) %*% as.matrix(data_gamma_freq)
        }
        
        gamma_a[gamma_a > nT] = nT
        gamma_v = data_gamma_freq/gamma_a
        
        
        ##
        dij = FDdistM
        dij = dij[rownames(dij) %in% rownames(data_2D), colnames(dij) %in% rownames(data_2D)]
        order_sp <- match(rownames(data_2D), rownames(dij))
        dij <- dij[order_sp, order_sp]
        
        dij = dij[rowSums(data_2D) > 0, rowSums(data_2D) > 0]
        data_2D = data_2D[rowSums(data_2D) > 0,]
        
        
        if (FDtau == 0) {
          dij[dij > 0] <- 1
          alpha_a = (1 - dij/1) %*% as.matrix(data_2D)
        } else {
          dij[which(dij > FDtau, arr.ind = T)] = FDtau
          alpha_a = (1 - dij/FDtau) %*% as.matrix(data_2D)
        }
        
        alpha_a[alpha_a > nT] = nT
        alpha_a = as.vector(alpha_a)
        alpha_v = rep(gamma_v, length(x))
        
        
        Chat = sapply(list(data_gamma, do.call(rbind, x)), function(y) iNEXT.3D:::Coverage(y, "incidence_raw", ncol(y)))
        mul_C2T = sapply(list(data_gamma, do.call(rbind, x)), function(y) iNEXT.3D:::Coverage(y, "incidence_raw", 2*ncol(y)))
        
        multiple = tibble('Assemblage' = c("Pooled assemblage", "Joint assemblage"), 
                          'T' = rep(ncol(data_gamma), 2), 
                          'U' = c(sum(data_gamma_freq), sum(data_2D)),
                          'S.obs' = c(sum(rowSums(data_gamma) > 0), sum(data_2D > 0)), 
                          'SC(T)' = Chat, 'SC(2T)' = mul_C2T,
                          'a1*' = c(sum(round(gamma_a) == 1), sum(round(alpha_a) == 1)), 
                          'a2*' = c(sum(round(gamma_a) == 2), sum(round(alpha_a) == 2)), 
                          'h1' = c(sum(gamma_v[round(gamma_a) == 1]), sum(alpha_v[round(alpha_a) == 1])), 
                          'h2' = c(sum(gamma_v[round(gamma_a) == 2]), sum(alpha_v[round(alpha_a) == 2])),
                          'Tau' = FDtau)
        
        
        single = DataInfo3D(x, diversity = "FD", datatype = "incidence_raw", FDtype = "tau_values", FDtau = FDtau, FDdistM = FDdistM)
        
        return(rbind(single, multiple) %>% cbind(Dataset = names(data)[i],.))
        
      }) %>% do.call(rbind,.)
      
    }
    
    rownames(output) = NULL
  }
  
  if (diversity == "FD" & FDtype == "AUC") {
    
    FDdistM = as.matrix(FDdistM)
    
    # if (datatype == 'abundance') {
    #   
    #   tmp <- lapply(data, rowSums)
    #   tmp = lapply(tmp, function(i) data.frame('value' = i) %>% rownames_to_column(var = "Species"))
    #   pdata = tmp[[1]]
    #   for(i in 2:length(tmp)){
    #     pdata = full_join(pdata, tmp[[i]], by = "Species")
    #   }
    #   pdata[is.na(pdata)] = 0
    #   pdata = pdata %>% column_to_rownames("Species")
    #   pdata = rowSums(pdata)
    #   
    #   order_sp <- match(names(pdata),rownames(FDdistM))
    #   FDdistM <- FDdistM[order_sp,order_sp]
    #   pdata <- matrix(pdata/sum(pdata), ncol = 1)
    #   
    # } else if (datatype == 'incidence_raw') {
    #   
    #   tmp <- lapply(data, function(x) {tmp = Reduce('+', x); tmp[tmp > 1] = 1; rowSums(tmp) })
    #   tmp = lapply(tmp, function(i) data.frame('value' = i) %>% rownames_to_column(var = "Species"))
    #   pdata = tmp[[1]]
    #   for(i in 2:length(tmp)){
    #     pdata = full_join(pdata, tmp[[i]], by = "Species")
    #   }
    #   pdata[is.na(pdata)] = 0
    #   pdata = pdata %>% column_to_rownames("Species")
    #   pdata = rowSums(pdata)
    #   
    #   order_sp <- match(names(pdata),rownames(FDdistM))
    #   FDdistM <- FDdistM[order_sp,order_sp]
    #   pdata <- matrix(pdata/sum(pdata), ncol = 1)
    #   
    # }
    # 
    # dmin <- min(FDdistM[FDdistM > 0])
    # dmean <- sum ( (pdata %*% t(pdata) ) * FDdistM) # dmean
    # dmax <- max(FDdistM[FDdistM > 0])
    # Tau = c(dmin, dmean, dmax)
    
    output = lapply(1:length(data), function(i) {
      
      x = data[[i]]
      
      if (datatype == 'abundance') {
        
        pdata = data.frame(rowSums(x))
        
        order_sp <- match(rownames(pdata),rownames(FDdistM))
        FDdistM <- FDdistM[order_sp,order_sp]
        pdata <- pdata/sum(pdata)
        
      } else if (datatype == 'incidence_raw') {
        
        tmp = Reduce('+', x)
        tmp[tmp > 1] = 1
        pdata = data.frame(rowSums(tmp))
        
        order_sp <- match(rownames(pdata),rownames(FDdistM))
        FDdistM <- FDdistM[order_sp,order_sp]
        pdata <- pdata/sum(pdata)
        
      }
      
      dmean <- sum ( (as.matrix(pdata) %*% t(as.matrix(pdata)) ) * FDdistM) # dmean
      
      dmin <- min(FDdistM[lower.tri(FDdistM)])
      dmax <- max(FDdistM[FDdistM > 0])
      
      if (datatype == "abundance") {
        
        if (is.null(colnames(x))) colnames(x) = paste('Assemblage_', 1:ncol(x), sep = "")
        
        mul_C2n = sapply(list("Pooled assemblage" = rowSums(x),
                              "Joint assemblage" = as.vector(x)), function(y) iNEXT.3D:::Coverage(y, "abundance", 2*sum(y)))
        
        multiple <- cbind(iNEXT.3D:::TDinfo(list("Pooled assemblage" = rowSums(x),
                                                 "Joint assemblage" = as.vector(x)), "abundance")[,1:5], 
                          "dmin" = dmin,
                          "dmean" = dmean,
                          "dmax" = dmax)
        
        single = DataInfo3D(x, diversity = "FD", datatype = "abundance", FDtype = "AUC", FDdistM = FDdistM)
        
        out = rbind(single, multiple) %>% cbind(Dataset = names(data)[i],.)
      }
      
      if (datatype == "incidence_raw") {
        
        if (is.null(names(x))) names(x) = paste('Assemblage_', 1:length(x), sep = "")
        
        data_gamma = Reduce('+', x)
        data_gamma[data_gamma > 1] = 1
        
        data_alpha = do.call(rbind, x)
        
        mul_C2T = sapply(list(data_gamma, data_alpha), function(y) iNEXT.3D:::Coverage(y, "incidence_raw", 2*ncol(y)))
        
        multiple <- cbind(iNEXT.3D:::TDinfo(list("Pooled assemblage" = c(ncol(data_gamma), as.vector(rowSums(data_gamma))),
                                                 "Joint assemblage" = c(ncol(data_alpha), as.vector(rowSums(data_alpha)))), "incidence_freq")[,1:6],
                          "dmin" = dmin,
                          "dmean" = dmean,
                          "dmax" = dmax)
        
        single = DataInfo3D(x, diversity = "FD", datatype = "incidence_raw", FDtype = "AUC", FDdistM = FDdistM)
        
        out = rbind(single, multiple) %>% cbind(Dataset = names(data)[i],.)
      }
      
      return(out)
      
      }) %>% do.call(rbind,.)
    
    rownames(output) = NULL
    }
  
  return(output)
  
  }




#' Printing iNEXTbeta3D object
#' 
#' \code{print.iNEXTbeta3D}: Print method for objects inheriting from class \code{iNEXTbeta3D}
#' @param x an \code{iNEXTbeta3D} object computed by \code{\link{iNEXTbeta3D}}.
#' @param ... additional arguments.
#' @export
print.iNEXTbeta3D <- function(x, ...){
  
  res <- lapply(x, function(y) {
    
    tmp = lapply(1:length(y), function(i) {
      
      # index1 = which(y[[i]]$SC == min(y[[i]]$SC))
      # index2 = which( (y[[i]]$Order.q != 0 & y[[i]]$SC == max(y[[i]]$SC)) | (y[[i]]$Order.q == 0 & y[[i]]$SC == max(y[[i]][y[[i]]$Order.q == 0, "SC"])) )
      # index3 = which(y[[i]]$Method %in% c("Observed", "Observed_SC(n, alpha)"))
      # 
      # index = sort(unique(c(index1, index1 + 1, index2, index2 - 1, index3 - 1, index3, index3 + 1)))
      # if (length(index) > 0) out = y[[i]][index,] else out = y[[i]]
      
      out = y[[i]]
      out[,c(3,5,7,8,9)] = round(out[,c(3,5,7,8,9)], 4)
      out[,4] = round(out[,4], 1)
      lapply(unique(out$Order.q), function(q) rbind(matrix(c("", paste("Order q = ", q, sep = ""), rep("", ncol(out)-2)), nrow = 1) %>% 
                                                      set_colnames(colnames(out)), 
                                                    out %>% filter(Order.q == q))) %>% do.call(rbind,.)
    })
    
    names(tmp) = names(y)
    
    return(tmp)
    }
  )
  
  print(res)
  # cat("\n")
  # cat("NOTE: The above output only shows seven rows for each table; call any table in iNEXT.beta3D.object to view complete output.\n", sep = "")
  return(invisible())
}





