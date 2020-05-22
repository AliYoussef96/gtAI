Introduction
=============

gtAI is a python package implemented to measure the tRNA adaptation index (tAI) [1], based on a novel approach [2]. The main advantages of the approach:
 
a) it is based only on the tRNA copy numbers (or tRNA levels) and the coding sequences, with no need for additional gene expression information.

b) it use a genetic algorithm, to optimize the sij-values-values (equation 1).

c) it overperform the two methods suggested before (the original tAI [1] and stAI [3]) to calculate the tAI.


**Note: the (g) in gtAI stand for genetic, as in the genetic algorithm, used in the implementation, but it calculates the know tAI [1].**