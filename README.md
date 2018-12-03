# DELMU
Simulation scripts used to produce the results presented in our paper R. LI et al. " DELMU: A Deep Learning Approach to Maximising the Utility of Virtualised Millimetre-Wave Backhauls"

This repository includes the 4 topologies we used to evaluate our proposed CNN approach to the non-convex optimisation problem, for each topology there are 2 folders containing
1) description of the topology in file capacity.npy and path.npy
2) matlab script to generate the global optimal solution via global search (GS), i.e. glbSearch.m
3) python script to solve the problem using greedy.py

For the purpose of reproducebility of our results, this repository also include the ramdom generated data trace for max and min flow demand in demandMaxPool.mat and demandMinPool.mat respectively.


