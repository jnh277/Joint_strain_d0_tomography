# Joint tomography of strain and the strain-free lattice spacing 
This repository contains Matlab code for two example simulations of jointly reconstruction the strain-free lattice spacing and the strain within a polycrystalline material from Bragg-edge neutron transmission measurements. The examples are a companion to the paper at 

The two examples discussed in the paper are included. These are;
> 'cantilever_beam_example.m': Reconstruction of the theoretical Saint-Venant strain field for a Cantilever beam
> 'C_shape_example.m': Reconstruction of a Finite Element Analysis (FEA) 'C' Shape sample subject to a compressive load
The examples can be found under the 'examples' folder.

The 'C' shape example uses pregenerated measurements. New sets of measurements can be generated using 'C_shape_simulate_measurements.m'

The examples are set up to run with hyper-parameters for the squared-exponential covariance function that were found in a previous optimisation. However, code to perform this optimisation is also included and can be enabled by setting the 'run_optimisation' parameter at the top of each script to 'true'. The optimisation procedure is very time consuming, it was found to take at least 10 hours on the workstation available to the authors.

It should be noted that the code has not been optimised for memory usage across different platforms. For instance on a 2018 macbook pro the cantilever beam example will run quite quickly without issues, whereas the 'C' shape example takes half an hour to run and uses all available memory.
