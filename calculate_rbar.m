%% Calculating average r, and its intra- and inter-chromosomal components

% This script uses the chromosome lengths (in chromosome_lengths.m) and crossover positions (in crossover_data.m)
% to calculate the average value of rbar, as well as its intra-chromosomal and inter-chromosomal components
% (from crossovers and independent assortment of homologs, respectively).

clc
clear

run crossover_data.m        % inputs the crossover position data (here, from Wang et al 2012)
run chromosome_lengths.m    % inputs chromosome lengths (here, as estimated from Fan et al 2011)

number_chroms = length(chrom_lengths);  % number of chromosomes in haploid set
relative_chrom_lengths = chrom_lengths/sum(chrom_lengths);  % Chromosomes' relative lengths sum to 1

number_cells = length(unique(data_mat(:,1)));   % number of cells from which the crossover data come

last_cell_label = max(data_mat(:,1));           % As in our data, the cells might not be labeled 1, ..., number_cells,
                                                % so number_cells does not equal last_cell_label

intra = zeros(number_cells,number_chroms);   % intra(m,n) will be the contribution of chromosome n to the
                        % intra-chromosomal component of average r in sperm m
                        % (there are 91 sperm cells, and 22 autosomes)
                        
crossovers = zeros(number_cells,1); % crossovers(m) will be the number of crossovers in sperm m

s = 0;  % s will enumerate each sequenced sperm, 1,2,...,91

identify_cell = zeros(number_cells,2); % identify_cell(s,j) links the enumeration of a sperm s with 
                             % its cell number in Wang et al 2012
identify_cell(:,1) = [1:number_cells]';

for i = 1:last_cell_label   % the cells in our data are enumerated from 1 to 144, but represent only 91 cells, 
                            % so we need to relabel them (according to s) from 1 to 91, in order
    
    if sum(data_mat(:,1)==i)==0 % If cell number i is not in the data set,
        s = s;                  % it doesn't count.
        
    elseif sum(data_mat(:,1)==i)>0  % If cell number i is in the data set,
        s = s+1;                    % it counts, and we proceed to calculate 
                                    % its average r.
                                    
        crossover_number = 0;   % crossover_number will sequentially count crossovers in this cell
                                    
        identify_cell(s,2) = i; 
        
        for c = 1:22
            
            if sum(data_mat(find(data_mat(:,1)==i),2)==c)==0 % If there are no crossovers on chromosome c,
                intra(s,c) = 0;                              % its contribution to intra-chromosomal recomb. is zero
                
            elseif sum(data_mat(find(data_mat(:,1)==i),2)==c)>0 % If there is at least one crossover along chrom. c,
                                                                % we proceed to calculate its contribution to the 
                                                                % intra-chromosomal component of average r
                                                                
                c_index = intersect(find(data_mat(:,1)==i), find(data_mat(:,2)==c));   % c_index is the set of indices of the 
                                                                                       % rows in the data matrix from file
                                                                                       % wang_et_al_2012_tableS2.m that 
                                                                                       % correspond to recombination events
                                                                                       % on chromosome c in the present sperm
                                                                                    
                crossover_number = crossover_number + length(c_index);  % adds the number of crossovers on the present chromosome
                                                                        % to our running tally of the number of crossovers in 
                                                                        % the present cell
                
                x_points = (data_mat(c_index,3) + data_mat(c_index,4))/2;   % we take the crossover points to be the
                                                                            % midpoints of the "start" and "end" sites of
                                                                            % crossover identification in the data
                                                                            
                chrom_with_crossovers = [0,x_points',chrom_lengths(c)];     % the start and end points of the chromosome,
                                                                            % and crossover positions in between (all in bp)
                                                                                           
                chrom_with_crossovers_rel = chrom_with_crossovers/chrom_lengths(c); % change positions such that this chroms'
                                                                                    % length normalized to 1
                
                odd = 0;    % sums lengths of haplotypes of one grandparentage
                even = 0;   % sums lengths of haplotypes of other grandparentage
                for marker = 1:length(chrom_with_crossovers)-1
                    if floor(marker/2)==marker/2
                        even = even + chrom_with_crossovers_rel(marker+1) - chrom_with_crossovers_rel(marker);
                    else
                        odd = odd + chrom_with_crossovers_rel(marker+1) - chrom_with_crossovers_rel(marker);
                    end
                end
                
                intra(s,c) = 2*even*odd*relative_chrom_lengths(c)^2;    % contribution to intra-chromosomal average r
                
            end
            
        end
        
        crossovers(s) = crossover_number;   % number of crossovers in cell s
        
    end
    
end

intras = sum(intra')';  % intra-chromosomal rbar for each cell

average_intra = sum(intras)/length(intras)  % average of the intra-chromosomal rbar across cells

stdev_intra = std(intras)   % standard deviation of intra-chromosomal rbar across cells

CoV_intra = stdev_intra/average_intra   % coefficient of variation

inter = 0.5*(1-relative_chrom_lengths'*relative_chrom_lengths) % inter-chrom. rbar, determined from chromosome lengths

rbar = average_intra + inter
                    
        