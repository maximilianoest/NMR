clc
clear all %#ok<CLALL>
close all
%% load data
configuration = readConfigurationFile('config.conf');
addpath(genpath(configuration.path2Library));
compartment = configuration.compartment;

disp(['Loading ' compartment ' data']);
data = loadResultsFromR1Simulation(configuration);

r1Perturbation = data.r1WithPerturbationTheory;
orientationAngles = rad2deg(data.orientationAngles);
positionAngles = rad2deg(data.positionAngles);
orientationsCount = length(orientationAngles);
positionsCount = length(positionAngles);

B0 = data.B0;
try
    startDateOfSimulation = data.startDateOfSimulation;
catch
    splittedFileName = split(fileName,'_');
    startDateOfSimulation = num2str(splittedFileName{1});
    msgbox(['Date: ' startDateOfSimulation],'Success');
end

%% set up system
B0WithoutComma = manipulateB0ValueForSavingPath(B0);
path2SaveFigures = [configuration.path2Results startDateOfSimulation ...
    '_RelaxationRatesOf' compartment 'At' ...
    num2str(B0WithoutComma) 'T.fig'];
path2SaveData = [configuration.path2Results startDateOfSimulation ...
    '_FilteredSimulationDataFrom' compartment 'At' ...
    num2str(B0WithoutComma) 'T.mat'];

try
    atomCount = data.atomCounter;
catch
    atomCount = findNumberOfCalculatedR1Rates(r1Perturbation);
end
calculatedR1Rates = r1Perturbation(:,:,1:atomCount);

%% trimm data
lowerPercentile = configuration.lowerPercentile;
upperPercentile = configuration.upperPercentile;
trimmedR1Rates = trimmR1DataWithPercentile(calculatedR1Rates, ...
    lowerPercentile,upperPercentile,orientationsCount,positionsCount);

atomCount = size(trimmedR1Rates,3);

%% calculate average and effective relaxation rates
meanRelaxationRates = mean(trimmedR1Rates,3);
medianRelaxationRates = median(trimmedR1Rates,3);

effectiveRelaxationRatesMean = mean(meanRelaxationRates,2);
% tmp = reshape(trimmedR1Rates,orientationsCount,atomCount*positionsCount);
% effectiveRelaxationRatesMedian = median(tmp,2);
effectiveRelaxationRatesMedian = mean(medianRelaxationRates,2);

%% saving data
data.calculatedR1Rates = calculatedR1Rates;
data.trimmedR1Rates = trimmedR1Rates;
data.upperPercentile = upperPercentile;
data.lowerPercentile = lowerPercentile;
data.meanRelaxationRates = meanRelaxationRates;
data.medianRelaxationRates = medianRelaxationRates;
data.effectiveRelaxationRatesMean = effectiveRelaxationRatesMean;
data.effectiveRelaxationRatesMedian = effectiveRelaxationRatesMedian;

save(path2SaveData,'-struct','data');

%% plotting data and saving them
fontSize = 14;
figs(1) = figure('DefaultAxesFontSize',fontSize);
plot(orientationAngles,effectiveRelaxationRatesMean,'LineWidth',1.5);
hold on
plot(orientationAngles,effectiveRelaxationRatesMedian,'LineWidth',1.5);
hold off
legend('Mean','Median','Location','NorthWest')
grid minor
title(['Dependency of Relaxation Rate from \theta at ' ...
    num2str(B0) ' Tesla (' num2str(lowerPercentile) '-' ...
    num2str(upperPercentile) '-Percentile)'])
xlabel('Angle \theta [°]')
ylabel('Relaxation Rates [Hz]')

figs(2) = figure('DefaultAxesFontSize',fontSize);
hold on
legendEntries = {};
legendEntryCounter = 1;
for positionNumber = 1:positionsCount
    positionAngle = positionAngles(positionNumber);
    plot(orientationAngles ...
        ,meanRelaxationRates(:,positionNumber),'LineWidth',1.5);
    legendEntries{legendEntryCounter} = ['Mean \phi: ' ...
        num2str(positionAngle)];  %#ok<SAGROW>
    legendEntryCounter = legendEntryCounter+1;
    
    plot(orientationAngles ...
        ,medianRelaxationRates(:,positionNumber),'LineWidth',1.5);
    legendEntries{legendEntryCounter} = ['Median \phi: ' ...
        num2str(positionAngle)];  %#ok<SAGROW>
    legendEntryCounter = legendEntryCounter+1;
end
hold off
legend(legendEntries,'Location','NorthWest')
title(['Relaxation Rates of ' compartment ])
xlabel('Angle \theta [°]')
ylabel('Relaxation Rates [Hz]')
grid minor

savefig(figs,path2SaveFigures)



