clc
clear all
close all


configuration = readConfigurationFile('config.conf');
addpath(genpath(configuration.path2Library));
compartment = configuration.compartment;

disp(['Loading ' compartment ' data']);
data = loadResultsFromR1Simulation(configuration);

r1Perturbation = data.r1WithPerturbationTheory;
atomCount = data.atomCounter;
calculatedR1Rates = r1Perturbation(2,:,1:atomCount);

lowerPercentile = configuration.lowerPercentile;
upperPercentile = configuration.upperPercentile;
trimmedR1Rates = trimmR1DataWithPercentile(calculatedR1Rates, ...
    lowerPercentile,upperPercentile,1,3);

meanR1 = mean(mean(trimmedR1Rates,3));
medianR1 = mean(median(trimmedR1Rates,3));



