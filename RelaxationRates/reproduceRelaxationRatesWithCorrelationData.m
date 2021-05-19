clc
clear all
%% load data
configuration = readConfigurationFile('config.conf');
addpath(genpath(configuration.path2Library));
path2ConstantsFile = configuration.path2ConstantsFileOnLocalMachine;

data = loadResultsFromR1Simulation(configuration);
B0 = data.B0;
startDateOfSimulation = data.startDateOfSimulation;
path2SaveFigures = [configuration.path2Results startDateOfSimulation ...
    '_ReconstructedRelaxationRates' configuration.compartment 'At' ...
    num2str(B0) 'T.fig'];

constants = readConstantsFile(path2ConstantsFile);
gammaRad = constants.gammaRad;
omega0 = gammaRad*B0;
hbar = constants.hbar;
Nm = constants.Nm;
mu0 = constants.mu0;
DD = 3/4*(mu0/(4*pi)*hbar*gammaRad^2)^2/(Nm^6);

[correlationFunctions1W0,correlationFunctions2W0] = ...
    getCorrelationFunctionsFromSimulationData(data);

atomCount = data.atomCounter;
deltaT = data.deltaT*data.shiftForCorrelationFunction;
lags = size(correlationFunctions1W0,4);

orientationAngles = rad2deg(data.orientationAngles);
positionAngles = rad2deg(data.positionAngles);
orientationCount = length(orientationAngles);
positionCount = length(positionAngles);

spectralDensities1W0 = zeros(orientationCount,positionCount,atomCount);
spectralDensities2W0 = zeros(orientationCount,positionCount,atomCount);
r1WithPerturbationTheory = zeros(orientationCount,positionCount);

for atomNumber = 1:atomCount
    for orientationNumber = 1:orientationCount
        for positionNumber = 1:positionCount
            [spectralDensities1W0(orientationNumber,positionNumber ...
                ,atomNumber),spectralDensities2W0(orientationNumber ...
                ,positionNumber,atomNumber)] = ...
                calculateSpectralDensities( ...
                squeeze(correlationFunctions1W0(orientationNumber ...
                ,positionNumber,atomNumber,:))' ...
                ,squeeze(correlationFunctions2W0(orientationNumber ...
                ,positionNumber,atomNumber,:))' ...
                ,omega0,deltaT,lags);
            
            r1WithPerturbationTheory(orientationNumber,positionNumber ...
                ,atomNumber) = calculateR1WithSpectralDensity( ...
                spectralDensities1W0(orientationNumber,positionNumber, ...
                atomNumber),spectralDensities2W0(orientationNumber ...
                ,positionNumber,atomNumber),DD);
        end
    end
end

%% plots
r1RatesMedian = median(r1WithPerturbationTheory,3);
r1RatesMean = mean(r1WithPerturbationTheory,3);

figs(1) = figure(1);
hold on
legendEntries = {};
legendEntryCounter = 1;
for positionNumber = 1:positionCount
    plot(orientationAngles,squeeze(r1RatesMean(:,positionNumber)) ...
        ,'LineWidth',1.5);
    legendEntries{legendEntryCounter} = ...
        ['Mean \phi: ' num2str(positionAngles(positionNumber)) '°']; %#ok<*SAGROW>
    legendEntryCounter = legendEntryCounter+1;
    
    plot(orientationAngles,squeeze(r1RatesMedian(:,positionNumber)) ...
        ,'LineWidth',1.5);
    legendEntries{legendEntryCounter} = ...
        ['Median \phi: ' num2str(positionAngles(positionNumber)) '°'];
    legendEntryCounter = legendEntryCounter+1;
end
grid minor
xlabel('Angle \theta [°]')
ylabel('Relaxation Rate [Hz]')
title('Relaxation Rates reproduced from Correlation Functions')
legend(legendEntries,'Location','NorthWest')
hold off

savefig(figs,path2SaveFigures);








