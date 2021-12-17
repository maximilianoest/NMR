clc
clear all
close all

configuration = readConfigurationFile('config.txt');
if configuration.runOnServer
    path2Results = configuration.path2ResultsOnServer;
    addpath(genpath(configuration.path2LibraryOnServer));
else
    path2Results = configuration.path2ResultsOnLocalMachine;
    path2Results = [path2Results 'nearestNeighboursAnalysis' '\']; 
    addpath(genpath(configuration.path2LibraryOnLocalMachine));
end

results = load("C:\Users\maxoe\Google Drive\Promotion\Results\nearestNeighboursAnalysis\Server\20211013_Results_relevantNearestNeighbours_water_H_50ns_05ps_wh.mat");

relaxationRates = results.r1WithPerturbationTheory; 
atomCounter = results.atomCounter;

averageR1 = squeeze(mean(relaxationRates(:,:,1:atomCounter,:),3));
effectiveR1 = squeeze(mean(averageR1,2));
experimentalData = squeeze(mean(effectiveR1,1));

medianR1 = squeeze(median(relaxationRates(:,:,1:atomCounter,:),3));
effectiveMedianR1 = squeeze(mean(medianR1,2));
overallMedianR1 = squeeze(mean(effectiveMedianR1,1));

try
    nearestNeighbourCases = results.nearestNeighbourCases;
catch
    nearestNeighbourCases = getValuesFromStringEnumeration( ...
        results.configuration.nearestNeighbourCases,';','numeric');
end
nearestNeighbourCases = flip(nearestNeighbourCases);

experimentalData = overallR1;
functionToFit = @(fittingParam,nearestNeighbourCases) ...
    log(nthroot(nearestNeighbourCases,3));
startingValues = [0];
fittedParameters = lsqcurvefit(functionToFit,startingValues ...
    ,nearestNeighbourCases,experimentalData);

figure(1)
hold on
% plot(experimentalData)
plot(log(nthroot(nearestNeighbourCases,3)))
hold off








