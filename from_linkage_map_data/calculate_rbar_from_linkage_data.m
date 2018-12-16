%% This script inputs chromosome lengths (chromosome_lengths.m) and linkage map data (via chromosome_indices.m) and 
%% then calculates, for a particular chromosome, the average contribution to the intra-chromosomal component of rbar.
%% Running the script for all chromosomes and summing the contributions will give the average value of the 
%% total intra-chromosomal component of rbar. 

clc
clear

female = 1; % 0 if male linkage map data, 1 if female

run chromosome_lengths

if female==1
load female_linkage_data.mat   % This loads an n x 4 matrix of inter-SNP map length data where, in row i, the 1st entry is 
                        % the chromosome number corresponding to a SNP, the 2nd entry is the SNP name, the 3rd entry is 
                        % the SNP position (HRG build 36 coordinates), the 4th entry is the estimated map distance (cM) 
                        % between the SNP in row i and that in row i-1 (if same chromosome).
elseif female == 0                                            
load male_linkage_data.mat
end

loci = 10^4;    % Number of evenly spaced pseudo-loci

chrom = 21;     % Chromosome whose contribution to the intra-chrom. value of rbar this script will calculate
               % Autosomes are 1-22; X chromosome is 23 (only for female data)
inds_chrom = find(linkage_data(:,1)==chrom); 
% Finds the row indices in the linkage map data matrix (linkage_data) that correspond to this chromosome

this_chrom_length = max(chrom_lengths(chrom), max(linkage_data(inds_chrom,3))+1);   % Chooses, as the chromosome length,
% the maximum of (i) the HRG chromosome length and (ii) the position of the last SNP in the linkage map data (in case 
% the latter is larger than the former). 
chrom_linkage_data = linkage_data(inds_chrom,:);    % Reduces the data to just the chromosome of interest
SNPs = chrom_linkage_data(:,3);                     % Just the SNP positions on this chromosome
distances = chrom_linkage_data(:,4)/100;            % Converts map distance scale from centiMorgans to Morgans

markers = round(linspace(0,this_chrom_length,loci));% Defines the positions of the evenly-spaced pseudomarkers

SNPs = [0;SNPs;this_chrom_length];                  % Adds 0 and the chrom length (or final SNP +1bp) to the SNP array

distances = [0;distances;0];                        % The average distances between every pair of SNPs, 
                                                    % including 0 and the end of the chromosome
cum_distances = cumsum(distances);                  % Cumulative distance as we move from the beginning to 
                                                    % the end of the chromosome, SNP by SNP

cum_distance_grid = zeros(loci,1);  cum_distance_grid(1) = 0;
below_inds = zeros(loci,1);         below_inds(1) = 1;
above_inds = zeros(loci,1);         above_inds(1) = 2;  % above_inds(i) and below_inds(i) will keep track, for each 
                                                        % pseudo-locus i, which SNP in the data is immediately before,
                                                        % and immediately above the pseudo-locus, moving from i = 1
                                                        % to i = loci. This will be stored as the index numbers of 
                                                        % these SNPs. 

lowerboundSNP = 1;  % Will keep track of which SNPs was immediately below the previous pseudo-locus,
                    % so that we don't have to look, for each pseudo-locus, at all SNPs to find out which SNP is
                    % immediately below the pseudo-locus (i.e., as we move from pseudo-locus to pseudo-locus, the 
                    % "lowerboundSNP" moves upwards. 
for marker = 2:loci
    distances =  SNPs([lowerboundSNP:end]) - markers(marker); % Distances (bp) between current pseudo-locus and all SNPs 
    above_ind = min(find(distances>-0.1))+lowerboundSNP-1;    % Which SNP is immediately above the current pseudo-locus?
    below_ind = above_ind - 1;                                % And which is immediately below?
    
    above_inds(marker) = above_ind; below_inds(marker) = below_ind; lowerboundSNP = below_ind;  % Reset lower bound SNP
                                                                                                % for next pseudo-locus
    
    above_ratio(marker) = (markers(marker) - SNPs(below_ind))/(SNPs(above_ind) - SNPs(below_ind));
    below_ratio(marker) = 1-above_ratio(marker);    % How far is the pseudomarker between the SNP immediately below  
                                                    % and immediately above it?
     
    cum_distance_grid(marker) =   below_ratio(marker)*cum_distances(below_ind)...
                                + above_ratio(marker)*cum_distances(above_ind); 
                                % The cumulative map distance up till the current pseudo-locus, calculated by
                                % interpolation from the cumulative map distance up till the SNP immediately below the
                                % pseudo-locus and the SNP immediately above it. 
    
end
    

intra_sum = 0;  % Will track the running total of the inra-chromosomal contribution to rbar as we move from locus
                % pair to locus pair
count_pairs = 0;% Will keep track of how many pairs we've counted---will eventually be loci x (loci-1)/2
for locus1 = 1:loci-1
    for locus2 = locus1+1:loci
        
        intra_sum = intra_sum + 0.5*tanh(2*(cum_distance_grid(locus2) - cum_distance_grid(locus1))); % Applies Kosambi's 
        % map function to the distance between pseudo-locus indexed "locus 1" and the pseudo-locus indexed "locus 2"
        count_pairs = count_pairs+1;
        
    end
end

if female == 0    % if male data
    intra_contrib = (intra_sum/count_pairs)*(this_chrom_length/sum(chrom_lengths(1:22)))^2

elseif female == 1  % if female data
    
    if chrom ~= 23  % if the current chromosome is not the X
    intra_contrib_noX = (intra_sum/count_pairs)*(this_chrom_length/sum(chrom_lengths(1:22)))^2
    end
    intra_contrib_X = (intra_sum/count_pairs)*(this_chrom_length/sum(chrom_lengths))^2

end