gtAI Description
=================

gtAI based on the optimization of Sij-values weights (equation 1) and relative codon synonymous usage (RSCU) values (equation 4).

The workflow to calculate gtAI:


1) Obtain a reference set by taking 5% (or more) of coding sequences from a tested genome, with the lowest effective number of codons values  (ENc) (equation 5).

2) Codons RSCU measured for this set of genes. 

3) The genetic algorithm used to search for the Sij-values that maximize the correlation between equation (1) and equation (4).

4) The gtAI calculated based on the optimized Sij-values from step 3.