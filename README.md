# miRcurio
Scripts for predicting miRNA using sRNA and degradome data.

For analysis miRcurio requires the perl package Temp and GetOpt::Std
Also requires local install of Bowtie, RNAfold and Patman

This program has two key parts. The first takes sRNA, maps them to the genome and returns a list of the most abundant sRNA mapping to a loci and their potential precursors. Putative precursors are folded and folds are checked to confirm they can generate a viable hairpin. Following this, sRNA are mapped to precursors to confirm if an miRNA* is in the dataset. 

The second part of this program takes putative miRNA from part one, takes the seed region of these and maps them to a degradome. Following degradome mapping the degradome sequences are mapped to the genome and the whole sRNA is mapped to the genome. Where the sRNA can map to regions where they have already had a seed match to the degradome they are output as potential miRNA.

This program returns the miRNA, the degradome cleavage site and target site in the genome, the hairpin and its genome coordinates.
