gtAI Workflow
==============


1) A reference set is obtained by taking 5% (or more) of coding sequences from a tested genome with the lowest ENc values (equation 3). (Or insert a reference set of interest)

2) Then, RSCU values for the reference set are generated (equation 4).

3) The genetic algorithm will search for the Sij-values that maximizes the correlation between (equation 1) and (equation 4).

4) The gtAI weights (Wi) are calculated based on the optimized Sij-values from step 3.

5) The weights are then normalized (wi), and gtAI for a coding sequence (g) is estimated using the (equation 3).
