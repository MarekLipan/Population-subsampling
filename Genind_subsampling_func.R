library(pacman)
p_load(ade4,adegenet,pegas,adegraphics, hierfstat)
#setwd("C:/Users/marek/Desktop/Verca")

# load the data
typ.genind <- read.structure(file="MSA_sheet_final_2018_selection_PCA.txt.stru", n.ind=202, n.loc=12, onerowperind=F, col.lab=1, col.pop=2, row.marknames=1, NA.char="-9", ask=FALSE) # set correct N individuals and loci
#typ.genind@pop

# transform into proper format
primula.hfstat <- genind2hierfstat(typ.genind)

#pe.statistics <- basic.stats(primula.hfstat,diploid=TRUE,digits=4)
#print(pe.statistics)
#pe.statistics[["n.ind.samp"]]

compute_hf_stat <- function(hf_data, S = 100){
  ############
  # Function for computing genetic statistics using subsampling from the populations.
  #
  # INPUTS:
  #
  # hf_data ...   DataFrame, obtained from genind using genind2hierfstat()
  #
  # S ...         Int, how many subsamples for each population we create
  #
  # OUTPUT:
  #
  # stat_mat ...  List, contains 2 matrices: 1) means (averages across subsamples and loci) of statistics "Ho", "Hs", "Fis"
  #               for each population 2) associated standard deviations
  #
  ############
  
  # names of statistics to be computed
  stat_names = c("Ho", "Hs", "Fis")

  # initialize matrix containing resulting statistics
  stat_mat_mean = matrix(data = 0, nrow = 3, ncol = nlevels(hf_data$pop), dimnames = list(stat_names, levels(hf_data$pop)))
  stat_mat_sd = matrix(data = 0, nrow = 3, ncol = nlevels(hf_data$pop), dimnames = list(stat_names, levels(hf_data$pop)))

  # the size of the subsamples is determined by the number of inidividual in the smallest population
  m = min(table(hf_data$pop))
  
  # subsampling
  for (s in 1:S){
    
    # initialize indices for the dataset with new subsamples
    sub_ind = numeric()
      
    # find the indices for the dataset with new subsamples (subsample from each populations)
    for (p in 1:max(as.numeric(hf_data$pop))){
      # add subsample indices for the p-th population
      sub_ind = c(sub_ind, sample(which(hf_data$pop == p) , m, replace=F))  #replace _ if true, the samples are used with replacement
    }
    
    
    # dataset with new subsamples
    hf_data_sub = hf_data[sub_ind,]
    
    # for the given subsamples, compute the list of full statistics
    full_stat_list = basic.stats(hf_data_sub,diploid=TRUE,digits=4)
    
    # fill the matrix of resulting statistics
    for (i in 1:3){
      stat = full_stat_list[[stat_names[i]]]
      # replace NaNs in the computed statistics by zeros
      stat[is.nan(stat)] = 0
      stat_mat_mean[i,] = stat_mat_mean[i,] + colMeans(stat, na.rm=T) # average across loci for each population
      stat_mat_sd[i,] = stat_mat_sd[i,] + apply(stat, 2, sd, na.rm=T) # apply sd to every column (population)
    }
  }

  # average across subsamples
  stat_mat_mean = stat_mat_mean/S
  stat_mat_sd = stat_mat_sd/S
  
  # combine results into list
  stat_mat = list(mean = stat_mat_mean, standard_deviation = stat_mat_sd)
  
  return(stat_mat)
}

# Example:
set.seed(123)
exm_results = compute_hf_stat(primula.hfstat, 100) #number of subsamples
exm_results$mean
exm_results$standard_deviation
