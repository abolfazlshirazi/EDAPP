-----------------------------------------------------
Constrained Estimation of Distribution Algorithm
By: Abolfazl Shirazi (ashirazi@bcamath.org)
Date: 2020-04-02
-----------------------------------------------------

This repository contains the codes for the algorithms
and the obtained data for the research on the constrained
Estimation of Distribution Algorithm. 

------------
FOLDERS
------------
    Algorithms:
        cmaes.m                 -> COVARIANCE MATRIX ADAPTATION EVOLUTIONARY STRATEGY
        geneticalgorithm.m      -> GENETIC ALGORITHM
        mgd.m                   -> CONSTRAINED ESTIMATION OF DISTRIBUTION ALGORITHM (Mixture of Gaussian Distribution)
        pso.m                   -> PARTICLE SWARM OPTIMIZATION

    InitialPop
        Contains initial seeds for comparing the algorithms.
        Seeds are generated for the benchmark of 13 problems.

    Problems:
        The benchmark of 13 constrained continuous optimization problems

    Results:
        The obtained results for algorithms:
            alg_mgd_lin_det         -> MGD + Linear Deterministic Method
            alg_mgd_lin_sto         -> MGD + Linear Stochastic Method
            alg_mgd_bis_det         -> MGD + Bisection Deterministic Method
            alg_mgd_bis_sto         -> MGD + Bisection Stochastic Method
            alg_cmaes               -> CMAES + Resampling Technique
            alg_ga                  -> GA + Penalty Function
            alg_pso                 -> PSO + Penalty Function
    
----------------------------------------------------
