%% Chromosome lengths

% This file inputs the lengths, in bp, of all chromosomes.

% We use as our estimate of the length of this particular individual's chromosomes
% the position of the last reported marker from Fan et al (2011)*, as listed in their Supplementary Data Set 4.
% *Fan et al (2011). Whole-genome molecular haplotyping of single cells. Nature biotechnology, 29(1), 51-57.

% We use this estimate instead of published chromosome lengths from reference genomes (e.g., HRG & UCSC)
% because in several cases, the last marker position in this individual exceeds 
% the length of the chromosome in the reference genomes (and, more particular to our purposes,
% some crossovers in this individual occur at positions beyond the length of the chromosome
% in the reference genomes.

chrom_lengths = [247175395
    242716953
    199375965
    191223332
    180630759
    170758410
    158819071
    146272231
    140221760
    135322654
    134446315
    132287411
    114118572
    106355404
    100338274
    88697596
    78637594
    76114852
    63786873
    62385570
    46920410
    49571074]; % lengths in bp, chromosomes ordered 1-22