These MATLAB scripts input chromosome lengths (chromosome_lengths.m) and crossover position data from individual gametes (crossover_data.m) and, using these, calculate the average value across gametes of the average genome-wide rate of genetic shuffling, rbar. They also calculate rbar's inter-chromosomal component (the contribution due to independent assortment of homologs) and the average value of its intra-chromosomal component (the contribution of crossovers). The metric rbar is described in Veller, Kleckner & Nowak (2018). There, the calculations that these scripts perform are carried out for crossover position data from 91 human sperm (data from Wang et al. 2012).

To carry out the calculations, run the file "calculate_rbar.m"


Veller, C., Kleckner, N., and Nowak, M.A. (2018). A rigorous measure of genome-wide genetic shuffling that accounts for crossover positions and Mendel's second law. biorXiv. doi.org/10.1101/194837

Wang, J., Fan, H.C., Behr, B., & Quake, S.R. (2012). Genome-wide single-cell analysis of recombination activity and de novo mutation rates in human sperm. Cell, 150(2), 402-412.
