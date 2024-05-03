# Optimal-scheduling-strategy-of-VPP
This is source code for paper (DOI: 10.16182/j.issn1004731x.joss.23-0840) published in Journal of System Simulation
<h3>About Documents</h3>
<h4>m 文件</h4>

Mainfunction.m is the main function and the functions are included:
 - Call the ScenarioMethod subfunction to get the example scenario after scenario generation and reduction, including the time, scenery output, and user load actual data
 - Scheduling model (constraints and objective functions) construction and solution
 - Numerical conversion, output results and plotting

ScenarioMethod.m is a subfunction, and the function includes:
 - Call the ManhattanDistance subfunction to get the scene set information after scene generation and cut, including the Manhattan distance matrix before and after cut, the probability distance matrix, the scene set, and the scene probability corresponding to the scene set
 - Obtaining typical scene probability and scenery output information under cutback conditions
 - Output results and plots

ManhattanDistance.m is a subfunction, and the function includes:
 - Calls the min12 subfunction, which handles the scene to be cut.
 - Implementation of the Manhattan distance method

min12.m is a subfunction, and the function includes:
 - Used to solve for the two smallest values and coordinates of a column or row vector.

<h2>Caveat</h2>
Regarding the <b>meteorological data</b> and <b>user load</b> two xlsx files can be imported into the project data according to the comments, this project these two raw data is not open source
