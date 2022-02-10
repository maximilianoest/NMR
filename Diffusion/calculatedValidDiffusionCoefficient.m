clc
clear all

% 1. find out how long the simulation time have to be
%   1.1: determine the average R1 for the whole simulation time
%   1.2: determine R1 for shorter simulation times and average over man
%   atoms
%   1.3: compare both to find the best simulation time to determine the
%   location of an atom
% 2. calculate the diffusion coefficient
%   2.1 % D = x^2/t


addpath(genpath('../library'));
configuration = readConfigurationFile('config.txt');
[path2Data,path2Save,path2ConstantsFile,path2LogFile] = ...
    setUpSystemBasedOnMachine(configuration);
deleteLogFile(path2LogFile);

logMessage('System is set up.', path2LogFile);

R1 = calculateR1ForWholeSimulationTime();




