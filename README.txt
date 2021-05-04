The constants are given in constants.txt and the configuration is given in config.conf.
For the constants the following units and explanations are given:

hbar 			= [Js] 
gammaRad		= [rad/Ts] gyromagnetic relation
B0 				= [T] main magnetic field
mu0				= [N/A^2] Vacuum permeability 4*pi*1e-7
Nm				= [m] Nanometer (Trajectory format)

The config file contains information about paths and computations. The names should be self
explaining.

Take care for the path of the library. If necessary take a look at the
code you try to run if you got an error that some function can't be found.


Full Datasets (used in paper):
- Myelin: /home/fschyboll/Dissertation/Paper4/Gromacs_Final/prd/Traj_Extracted
	-> 250ns simulation time and time between sampling points 2ps
	-> for simulation the step size of this data set must be changed from 1 to 2 and thus from 1ps to 2 ps.
	-> example trajectory has already 250ns and 2ps step size (but it can't be used bacause it isn't representative)

- water: /home/fschyboll/Dissertation/Paper3/Matlab/Eval_MonoLayer3/Traj_mat/Short_50ns_025ps/water_H_50ns_025ps_wh.mat
	-> 50ns simulation time and time between sampling points 0.25ps


