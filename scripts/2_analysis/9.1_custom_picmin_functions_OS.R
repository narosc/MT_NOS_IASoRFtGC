### This function is based on the picmin vignette: https://github.com/TBooker/PicMin/tree/main/vignettes

### Some comments have been added for clarity, plus integration with parallel library for computational speed

## Load libraries
library(PicMin)
library(foreach)
library(doParallel)

###
##### This is how PicMin generates the null distribution of p-values for the number of species of interst (n) used in this iteration
###

PicMinNull = function(linMin, linMax) {


  nullP <- foreach(n=linMin:linMax) %dopar% {
    
  # Run 10,000 replicate simulations of this situation and build the correlation matrix
  emp_p_null_dat <- t(replicate(40000, PicMin:::GenerateNullData(1.0, n, 0.5, 3, 10000)))
  # Calculate the order statistics p-values for each simulation
  emp_p_null_dat_unscaled <- t(apply(emp_p_null_dat ,1, PicMin:::orderStatsPValues))
  # Use those p-values to construct the correlation matrix
  null_pMax_cor_unscaled <- cor( emp_p_null_dat_unscaled )
  
  return(null_pMax_cor_unscaled)
  
  }
  
  names(nullP) = paste0('l',linMin:linMax)
  
  return(nullP)
}


## Slight modification of the PicMin pipeline. Key differences: 
#1) the null P distribution is pre-calculated in an external function (above), so that it doesn't have to be re-calculated every time we run the script for the same number of species. 
#2) the number of repetition to calculate P values is set iteratively. First, a standard round with a relatively low number of repetitions (e.g. 1000), applied to all genes. This allows to estimate p with a rough accuracy above the number of permutations (P=1-0.001). If p are smaller than the number of permutations (P<0.001),  we iteratively re-calculate P by running x10 more permutations. This continues until a precise small p-value is found, or until reaching a maximal number of permutations allowed (e.g. 10^9), to prevent infinite loops. The advantage of this method: getting very precise small p-values, which are important to find significant SNPs after q-value correction.   

RunPicmin = function(all_lins_p, ### this a table where every row is a genomic window/gene, every column a species, and the values are the p-values of the geas
                     numReps_std=1000, ### number of repetitions applied to all genes to calculate rough p-value
                     numReps_max=10^9, ### maximal number of repetitions allowed to potentially significant genes to calculate exact p-value 
                     minObsData=3, ### this is a number indicating the minimal number of species used to search for convergence. The guidelines recommend >=3. Note: this is the number of species with non missing values, so the p-values table could have 10 species, but  only those with a least 3 species with data (vs. 7 with missing values) will be used. 
                     nullP) {  ### this is a list object storing the null distribution of p-values for different number of lineages

   ### Set run parameters
  nLins=ncol(all_lins_p) # calculate how many lineages/species are there in the dataset (how many columns in the p-value matrix)
  missingDataLevels=c(minObsData:nLins) # set all the missing data levels that will be checked for overlap. For ex. if this is 3:10, it will check overlap for genes that are non-missing in 3 species, then 4, 5, and so on until 10. 

  ### Run loop and store results for every missing data level
  results = list()
  count = 0 
  
  for (n in missingDataLevels){ # for every missing data level

   
    count = count + 1
    
    ###
    ##### Then we only retain genes that are observed in exactly the number of species of interest (n) in this iteration 
    ###
    # Screen out gene with no evidence for adaptation
    lins_p_n <-  as.matrix(all_lins_p[rowSums(is.na(all_lins_p)) == nLins-n,])
    
    
    ###
    ##### Finally, we run PicMin gene-by-gene. For every gene there are 2 outputs: an p-value of the rank analysis (compared to null distribution expectation), and the n_est parameter, indicating how many species have a p-value with ranks lower than expected by chance 
    ###
    
    if (dim(lins_p_n)[1] ==0){
      next
    }

    
    res <- foreach(i=seq(nrow(lins_p_n)), .combine=rbind) %dopar% {

      nri = numReps_std # number of permutation used for this gene
      
      ## try test with standard number of permutations
      test_result <- PicMin:::PicMin(na.omit(lins_p_n[i,]), nullP[[paste0('l',n)]], numReps = nri)
      res_p <- test_result$p
      res_n <- test_result$config_est
      
       ## if p-value is minimal for current number of permutations --> increase number of permutations....
      while (res_p < 1/nri & nri < numReps_max) { # ...unless maximal number of permutations is reached (to avoid infinite loops)
        
        nri = nri*10
        
        ## try test with increased number of permutations
        print('doing it')
        test_result <- PicMin:::PicMin(na.omit(lins_p_n[i,]), nullP[[paste0('l',n)]], numReps = nri)
        res_p <- test_result$p
        res_n <- test_result$config_est
        
      }
      
      return(data.frame(res_p, res_n))
    }
    
    
    results[[count]] = data.frame(numLin = n,
                                  p = res$res_p,
                                  q = p.adjust(res$res_p, method = "fdr"),
                                  n_est = res$res_n,
                                  locus = row.names(lins_p_n) )
    
  }
  
  ### Results are formatted in a table
  picMin_results <- do.call(rbind, results)
  
  ### P-values across runs are corrected using FDR
  picMin_results$pooled_q <- p.adjust(picMin_results$p, method = "fdr")

  return(picMin_results)
  
}




