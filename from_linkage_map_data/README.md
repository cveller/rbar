The MATLAB script loads a data file containing map lengths between adjacent SNPs along all human chromosomes, for either females ("female_linkage_data.mat") or males ("male_linkage_data.mat"), and calculates from these data the average intra-chromosomal contribution to rbar for a particular chromosome. Carrying out the calculation for all chromosomes, and summing the results, gives the total average intra-chromosomal contribution to rbar. The data used are from Kong et al. (2010). The metric rbar is described in Veller, Kleckner & Nowak (2018), which reports the results of the calculations performed by this script on the data from Kong et al. (2010).

To carry out the calculation, run the file "calculate_rbar_from_linkage_data.m", specifying in the script the chromosome number and whether to run the script for the male or female data.



Kong, A., Thorleifsson, G., Gudbjartsson, D.F., Masson, G., Sigurdsson, A., Jonasdottir, A., ... & Gudjonsson, S.A. (2010). Fine-scale recombination rate differences between sexes, populations and individuals. Nature, 467(7319), 1099-1103. 

Veller, C., Kleckner, N., and Nowak, M.A. (2018). A rigorous measure of genome-wide genetic shuffling that accounts for crossover positions and Mendel's second law. biorXiv. doi.org/10.1101/194837
