# EDA++
Estimation of Distribution Algorithms with Feasibility Conserving Mechanisms for Constrained Continuous Optimization

## Suplementary Materials
This repository contains the codes for the algorithms and the obtained data for the research on the constrained Estimation of Distribution Algorithm EDA++.

### Analysis of the components of the algorithm

* **Algorithms**

  - cmaes.m                 : COVARIANCE MATRIX ADAPTATION EVOLUTIONARY STRATEGY
  - geneticalgorithm.m      : GENETIC ALGORITHM
  - mgd.m                   : CONSTRAINED ESTIMATION OF DISTRIBUTION ALGORITHM (Mixture of Gaussian Distribution)
  - pso.m                   : PARTICLE SWARM OPTIMIZATION

* **InitialPop**

    Contains initial seeds for comparing the algorithms. Seeds are generated for the benchmark of 13 problems.

* **Problems**

    The benchmark of 13 constrained continuous optimization problems.

* **Results**

    The obtained results for algorithms:
    - alg_mgd_lin_det         : MGD + Linear Deterministic Method
    - alg_mgd_lin_sto         : MGD + Linear Stochastic Method
    - alg_mgd_bis_det         : MGD + Bisection Deterministic Method
    - alg_mgd_bis_sto         : MGD + Bisection Stochastic Method
    - alg_cmaes               : CMAES + Resampling Technique
    - alg_ga                  : GA + Penalty Function
    - alg_pso                 : PSO + Penalty Function

### Comparison with other algorithms

* **Algorithms**

  - ALG_BPMAGES                 : BP-ϵMAg-ES 
  - ALG_COLSHADE                : COLSHADE
  - ALG_EDAPP                   : EDA++
  - ALG_ENMODE                  : EnMODE
  - ALG_LSHADE44                : LSHADE44

* **CEC2020**

    CEC 2020 real-world constrained optimization benchmark problems.

* **Solutions**

    The obtained solutions for each algorithm regarding all problems in the benchmark (5 algorithms, 57 problems, 25 runs).
  
## Credit

Abolfazl Shirazi

Basque Center for Applied Mathematics - BCAM | Universidad del Pais Vasco UPV/EHU

Email: ashirazi@bcamath.org

Web: http://www.homasim.com/ashirazi


