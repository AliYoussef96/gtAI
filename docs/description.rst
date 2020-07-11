gtAI Workflow
==============


1) A reference set is obtained by taking 5% (or more) of coding sequences from a tested genome with the lowest ENc values (equation 3). (Or insert a reference set of interest)

2) Then, RSCU values for the reference set are generated (equation 4).

3) The genetic algorithm will search for the Sij weights that maximizes the correlation between (equation 1) and (equation 4).

4) The final Wi values are calculated based on the optimized Sij weights from step 3. 

5) Calculated Wi values are normalized (wi), then the tAI value for a coding sequence (g) is calculated using equation 3.
