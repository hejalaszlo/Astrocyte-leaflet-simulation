# Simulation of Ca2+ dynamics in astrocyte leaflets
Matlab class and helper functions to simulate astrocytic calcium and sodium dynamics in thin astrocytic processes of genuine, experimentally obtained geometries.

To simulate Ca2+ dynamics in experimentally determined synapses, do the followings:

1)	Add the @Synapse folder to your path in Matlab

2)	Download our pre-processed synapse data based on the Kasthuri et al., 2015. It contains 1700 selected synapses.
http://downloadables.ttk.hu/heja/Front_Cell_Neurosci_2021/vastSynapses.mat

3)  Load this file to Matlab:
load("vastSynapses.mat");

4)	Now, we will run a simulation on the first synapse in this database.

5)	Since tissue is shrinked in the EM samples, we need to correct for the shrinkage:
synapses(1).correctECS;

6)	Clear previous simulation data:
synapses(1).simulationClear;

7)	Run the simulation:
synapses(1).simulationRun;

8) The simulation results can be viewed at synapses(1).Simulation or you can export it to a .mat file:
synapses(1).simulationExport("synapse1_simulation");

# Citation
If you find this tool useful, please, cite the following publication:
Heja L, Szabo Z, Peter M, Kardos J. Spontaneous Ca2+ Fluctuations Arise in Thin Astrocytic Processes With Real 3D Geometry. Front Cell Neurosci., 15:617989, 2021
https://doi.org/10.3389/fncel.2021.617989