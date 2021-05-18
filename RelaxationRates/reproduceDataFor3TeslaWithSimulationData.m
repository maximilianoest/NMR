clc 
clear all
%% load data
configuration = readConfigurationFile('config.conf');
addpath(genpath(configuration.path2Library));
path2ConstantsFile = configuration.path2ConstantsFileOnLocalMachine;

B0 = 3;
constants = readConstantsFile(path2ConstantsFile);
gammaRad = constants.gammaRad;
omega0 = gammaRad*B0;
hbar = constants.hbar;
Nm = constants.Nm;
mu0 = constants.mu0;
DD = 3/4*(mu0/(4*pi)*hbar*gammaRad^2)^2/(Nm^6);

data = loadResultsFromR1Simulation(configuration);

sumCorrelationFunction1W0 = data.sumCorrelationFunctionW0Saver;
sumCorrelationFunction2W0 = data.sumCorrelationFunction2W0Saver;
atomCount = data.atomCounter;
deltaT = data.deltaT;
lags = data.lags;

orientationCount = size(sumCorrelationFunction1W0,1);
positionCount = size(sumCorrelationFunction1W0,2);

correlationFunction1W0 = sumCorrelationFunction1W0/atomCount;
correlationFunction2W0 = sumCorrelationFunction2W0/atomCount;

spectralDensity1W0 = zeros(orientationCount,positionCount);
spectralDensity2W0 = zeros(orientationCount,positionCount);
r1WithPerturbationTheory = zeros(orientationCount,positionCount);

for orientationNumber = 1:orientationCount
    for positionNumber = 1:positionCount
        [spectralDensity1W0(orientationNumber,positionNumber) ...
            ,spectralDensity2W0(orientationNumber,positionNumber)] = ...
            calculateSpectralDensities( ...
            correlationFunction1W0(orientationNumber,positionNumber) ...
            ,correlationFunction2W0(orientationNumber,positionNumber) ...
            ,omega0,deltaT,lags);
        
        r1WithPerturbationTheory(orientationNumber,positionNumber) ...
            = calculateR1WithSpectralDensity( ...
            spectralDensity1W0(orientationNumber,positionNumber) ...
            ,spectralDensity2W0(orientationNumber,positionNumber),DD);
    end
end


%% plots

figure(1)







