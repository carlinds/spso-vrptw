# Set-based Particle Swarm Optimization for the VRPTW
This project is a Python-implementation of a Set-based Particle Swarm Optimization (SPSO) for solving the Vehicle Routing Problem with Time Windows (VRPTW). 
The SPSO is implemented purely for testing purposes, more specifically for testing on the [Solomon instances](https://www.sintef.no/projectweb/top/vrptw/solomon-benchmark/). The results from these tests are presented below and discussed in more detail in the paper "Comparison of Exact and Approximate methods for the Vehicle Routing Problem with Time Windows".

## Usage
A single Solomon-instance can be evaluated using `spso\evaluateInstance.py`. All tuning parameters and run-time options are specified in this file.
All instances from the Solomon benchmarks are found in `solomon-instances`.

## Results
The results used in the paper are stored under `results`, where `results\solomon-instances` are the results from evaluating all instances and `results\global-updates` are the results from evaluating a few selected instances with global updates enabled. A more visual representation of the results are presented below.

### Solomon benchmark
Results from the SPSO on the Solomon instances, compared to a MILP algorithm. For more details and discussion regarding the results, see the paper.

#### 25 customers
<img src="https://github.com/carlinds/spso-vrptw/blob/master/results/milp_comparison_25_customers.png" alt="Solomon benchmark for 25 customers" width="700"/>

#### 50 customers
<img src="https://github.com/carlinds/spso-vrptw/blob/master/results/milp_comparison_50_customers.png" alt="Solomon benchmark for 50 customers" width="700"/>

#### 100 customers
<img src="https://github.com/carlinds/spso-vrptw/blob/master/results/milp_comparison_100_customers.png" alt="Solomon benchmark for 100 customers" width="700"/>
