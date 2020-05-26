Introduction
=============

gtAI is a package implemented in python to measure the tRNA adaptation index (tAI) [1], based on a novel approach [2]. The main advantages of this approach:
 
1) It is based on the tRNA gene copy number (or tRNA levels) and the coding sequences, with no need for additional gene expression information. (can be used If available)

2) It uses a genetic algorithm to optimize the sij-values (equation 1).

3) It outperforms previously suggested methods (the original tAI [1] and stAI [3]) to calculate the tAI.


**Note: The "g" in gtAI stands for genetic, as in the genetic algorithm used in the implementation, but it calculates the known tAI [1].**
