clc
clear all %#ok<CLALL>
%% load data
configuration = readConfigurationFile('config.conf');
addpath(genpath(configuration.path2Library));

data2Load = [configuration.path2Data configuration.fileName '.mat'];
data = load(data2Load);

r1Perturbation = data.r1WithPerturbationTheory;
orientationAngles = rad2deg(data.orientationAngles);
positionAngles = rad2deg(data.positionAngles);
orientationsCount = length(orientationAngles);
positionsCount = length(positionAngles);

B0 = data.B0;
dateFromConfig = num2str(configuration.date);

%% set up system
B0WithoutComma = manipulateB0ValueForSavingPath(B0);
path2SaveFigures = [configuration.path2Results dateFromConfig ...
    'RelaxationRatesOf' configuration.compartment 'At' ...
    num2str(B0WithoutComma) 'T.fig'];
path2SaveData = [configuration.path2Results dateFromConfig ...
    '_ResultsFrom' configuration.compartment '.mat'];

try
    atomCount = data.atomCounter;
catch
    atomCount = findNumberOfCalculatedR1Rates(r1Perturbation);
end

calculatedR1Rates = r1Perturbation(:,:,1:atomCount);

%% trimm data
lowerPercentile = 30;
upperPercentile = 100-lowerPercentile;
trimmedR1Rates = trimmR1DataWithPercentile(calculatedR1Rates, ...
    lowerPercentile,upperPercentile,orientationsCount,positionsCount);

atomCount = size(trimmedR1Rates,3);

%% calculate average and effective relaxation rates
meanRelaxationRates = mean(trimmedR1Rates,3);
medianRelaxationRates = median(trimmedR1Rates,3);

effectiveRelaxationRatesMean = mean(meanRelaxationRates,2);
tmp = reshape(trimmedR1Rates,orientationsCount,atomCount*positionsCount);
effectiveRelaxationRatesMedian = median(tmp,2);

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
figs(1) = figure(1);
plot(orientationAngles,effectiveRelaxationRatesMean,'LineWidth',1.5);
hold on
plot(orientationAngles,effectiveRelaxationRatesMedian,'LineWidth',1.5);
hold off
legend('Mean','Median')
grid on
title(['Angle Dependency of Relaxation Rate at ' num2str(B0) ' Tesla (' ...
    num2str(lowerPercentile) '-'  num2str(upperPercentile) '-Percentile)'])
xlabel('Angle [�]')
ylabel('Relaxation Rates [Hz]')

savefig(figs,path2SaveFigures)


