Introduction
=============

gtAI is a package implemented in python to measure the tRNA adaptation index (tAI) [1], based on a novel approach [2]. The main advantages of this approach:
 
1) It requires the tRNA gene copy number (or tRNA levels) and coding sequences, without the need for additional gene expression information (can be used If available).

2) It uses a genetic algorithm to reach the best set of Sij-values (equation 1).

3) It outperforms previously suggested methods (the original tAI [1] and stAI [3]) in tRNA adaptation index (tAI) computation by producing significantly better results.


**Note: The "g" in gtAI stands for genetic, as in the genetic algorithm used in the implementation, but it calculates the known tAI [1].**
