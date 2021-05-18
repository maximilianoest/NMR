
%% load data
clc
clear all
close all

configuration = readConfigurationFile('config.conf');
addpath(genpath(configuration.path2Library));
compartment = configuration.compartment;

data = loadResultsFromR1Simulation(configuration);

atomCount = data.atomCounter;
orientationAngles = rad2deg(data.orientationAngles);
orientationsCount = length(orientationAngles);
positionAngles = rad2deg(data.positionAngles);
positionsCount = length(positionAngles);
deltaT = data.deltaT*data.shiftForCorrelationFunction;

B0 = data.B0;
B0WithoutComma = manipulateB0ValueForSavingPath(B0);

startDateOfSimulation = data.startDateOfSimulation;
path2SaveFigures = [configuration.path2Results startDateOfSimulation ...
    '_CorrelationFunctionsOf' compartment 'At' ...
    num2str(B0WithoutComma) 'T.fig'];
path2SaveData = [configuration.path2Results startDateOfSimulation ...
    '_ResultsFrom' compartment 'At' num2str(B0WithoutComma) 'T.mat'];

%% analyze data
correlationFunctions1W0 = data.correlationFunction1W0Saver;
correlationFunctions2W0 = data.correlationFunction2W0Saver;
atomIndex = data.atomIndex;
atomIndex = atomIndex(1:atomCount);
correlationFunctions1W0 = ...
    squeeze(correlationFunctions1W0(:,:,atomIndex,:));
correlationFunctions2W0 = ...
    squeeze(correlationFunctions2W0(:,:,atomIndex,:));

averageCorrelationFunctions1W0 = squeeze(mean( ...
    correlationFunctions1W0,3));
averageCorrelationFunctions2W0 = squeeze(mean( ...
    correlationFunctions2W0,3));

effectiveCorrelationFunctions1W0 = squeeze(mean( ...
    averageCorrelationFunctions1W0,2));
effectiveCorrelationFunctions2W0 = squeeze(mean( ...
    averageCorrelationFunctions2W0,2));

%% plotting
tauAxis = 0:deltaT:(size(correlationFunctions1W0,4)-1)*deltaT;
tauMin = 0;
tauMax = 5e-7;
valueMin = 0;
valueMax = 0.5;
fontSize = 16;

legendEntries = {};
figs(1) = figure('DefaultAxesFontSize',fontSize);
hold on
for orientationNumber = 1:orientationsCount
    plot(tauAxis,abs(real(effectiveCorrelationFunctions1W0( ...
        orientationNumber,:)/effectiveCorrelationFunctions1W0( ...
        orientationNumber,1))),'LineWidth',1.5)
    legendEntries{orientationNumber} = num2str( ...
        orientationAngles(orientationNumber),'Theta: %.2f'); %#ok<SAGROW>
    axis([tauMin tauMax valueMin valueMax])
end
hold off
grid minor
legend(legendEntries)
title('Effective Correlation Functions for w0')
xlabel('\tau')

legendEntries = {};
figs(2) = figure('DefaultAxesFontSize',fontSize);
hold on
for orientationNumber = 1:orientationsCount
    for positionNumber = 1:positionsCount
        plot(tauAxis,squeeze(abs(real(averageCorrelationFunctions1W0( ...
            orientationNumber,positionNumber,:) ...
            /averageCorrelationFunctions1W0( ...
            orientationNumber,positionNumber,1)))),'LineWidth',1.5)
    end
    legendEntries{orientationNumber} = num2str( ...
        orientationAngles(orientationNumber),'Theta: %.2f'); %#ok<SAGROW>
end
axis([tauMin tauMax valueMin valueMax])
hold off
grid minor
legend(legendEntries)
title('Average Correlation Functions for w0')
xlabel('\tau')

legendEntries = {};
figs(3) = figure('DefaultAxesFontSize',fontSize);
hold on
for orientationNumber = 1:orientationsCount
    plot(tauAxis,abs(real(effectiveCorrelationFunctions2W0( ...
        orientationNumber,:)/effectiveCorrelationFunctions2W0( ...
        orientationNumber,1))),'LineWidth',1.5)
    legendEntries{orientationNumber} = num2str( ...
        orientationAngles(orientationNumber),'Theta: %.2f'); %#ok<SAGROW>
    axis([tauMin tauMax valueMin valueMax])
end
hold off
grid minor
legend(legendEntries)
title('Effective Correlation Functions for 2w0')
xlabel('\tau')

legendEntries = {};
figs(4) = figure('DefaultAxesFontSize',fontSize);
hold on
for orientationNumber = 1:orientationsCount
    for positionNumber = 1:positionsCount
        plot(tauAxis,squeeze(abs(real(averageCorrelationFunctions2W0( ...
            orientationNumber,positionNumber,:) ...
            /averageCorrelationFunctions2W0( ...
            orientationNumber,positionNumber,1)))),'LineWidth',1.5)
    end
    legendEntries{orientationNumber} = num2str( ...
        orientationAngles(orientationNumber),'Theta: %.2f'); %#ok<SAGROW>
end
axis([tauMin tauMax valueMin valueMax])
hold off
grid minor
legend(legendEntries)
title('Average Correlation Functions for w0')
xlabel('\tau')
ylabel('G(\tau)')

savefig(figs,path2SaveFigures);




