The constants are given in constants.txt and the configuration is given in config.conf.
For the constants the following units and explanations are given:

hbar 			= [Js] 
gammaRad		= [rad/Ts] gyromagnetic relation
B0 			= [T] main magnetic field
mu0			= [N/A^2] Vacuum permeability 4*pi*1e-7
Nm			= [m] Nanometer (Trajectory format)
DD 			= [J/rad] Dipol-Dipol ineraction constants (Calculated in Script)
omega0			= [rad/s] Larmor (anglular) frequency (Calculated in Script)

The config file contains information about paths and computations. The names should be self
explaining.

Processing Stream:
1) calculate the R1 rates with the orientation dependency script
2) with the script analyzeAngleDependency.m calculate the trimmed values for the R1 rates and
	calculate the data for the prediction of the R1 rate
3) with the script MAIN_calculateR1RatesWithOptimizedParameters.m calculate the predicted R1 rate.



