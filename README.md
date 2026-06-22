# ZCZU for CEC 2026 CMOP

## Authors

* Changzhe Zheng

Contact: Z2377593163@163.com

## Algorithm Name

ZCZU

## Algorithm Explanation

ZCZU is a MATLAB-based constrained multi-objective evolutionary algorithm implemented on the PlatEMO platform. The algorithm uses a dual-population framework with dynamic constraint processing and resource allocation to solve constrained multi-objective optimization problems.

This repository provides the source code of ZCZU for the CEC 2026 CMOP competition.

## Pseudocode
Algorithm: ZCZU

Input:
    Problem, population size N, maximum function evaluations Max_FEs

Output:
    Final population

1. Initialize two populations, Population1 and Population2.
2. Evaluate the objective values and constraint violations of both populations.
3. Initialize the fitness values, archive, and constraint-processing status.

4. While the termination condition is not satisfied:
      4.1 Update the search-status information.
      4.2 Determine whether the algorithm should enter the next constraint-processing stage.
      4.3 Generate offspring for the two populations.
      4.4 Update the external archive.
      4.5 Perform environmental selection for Population1.
      4.6 Perform environmental selection for Population2.
      4.7 Update the resource allocation ratio according to offspring success rates.
      4.8 Record Min_IGD and MCV.

5. Return the final population.


## Experimental Settings

* Problems: SDC1-SDC15
* Runs: 30 independent runs
* Max_FEs: 200000
* Dimension: D = 30
* Population size: 100
* Indicators: Min_IGD and MCV
* Platform: MATLAB + PlatEMO
* Training: No

## Usage

Please place the ZCZU folder into the algorithm directory of PlatEMO and run the algorithm on SDC1-SDC15 using the settings above.

