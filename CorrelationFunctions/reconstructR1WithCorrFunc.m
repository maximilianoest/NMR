clc
clear all
close all

results = load("C:\Users\maxoe\Google Drive\Promotion\Results\relevantNearestNeighboursCorrelationFunction\20220422_Results_relevantNearestNeighboursCorrelationFunction_DOPSlipid4.mat");
baseConfiguration = readConfigurationFile('../baseConfiguration.txt');
addpath(genpath(baseConfiguration.path2LibraryOnLocalMachine));

orientationAngles = results.orientationAngles;
positionAngles = results.positionAngles;

disp('-- Only summed up correlation function --')
corrFunc1W0 = results.sumCorrelationFunction1W0Saver/results.atomCounter;
corrFunc2W0 = results.sumCorrelationFunction2W0Saver/results.atomCounter;

for orientationNr = 1:length(orientationAngles)
    for positionNr = 1:length(positionAngles)
        [spectralDensity1W0, spectralDensity2W0] = ...
            calculateSpectralDensities(squeeze(corrFunc1W0(orientationNr ...
            ,positionNr,:))',squeeze(corrFunc2W0(orientationNr,positionNr,:))' ...
            ,results.omega0,results.samplingFrequency,results.lags);
        r1WithPerturbationTheory(orientationNr,positionNr) = ...
            calculateR1WithSpectralDensity(spectralDensity1W0 ...
            ,spectralDensity2W0,results.dipolDipolConstant);
    end
end
r1WithPerturbationTheory
r1 = mean(results.r1WithPerturbationTheory(:,:,:,1),3)

corrFunc1W0 = results.correlationFunction1W0Saver;
corrFunc2W0 = results.correlationFunction2W0Saver;

disp('-- Averaging already in simulation --')

for orientationNr = 1:length(orientationAngles)
    for positionNr = 1:length(positionAngles)
        [spectralDensity1W0, spectralDensity2W0] = ...
            calculateSpectralDensities(squeeze(corrFunc1W0(orientationNr ...
            ,positionNr,:))',squeeze(corrFunc2W0(orientationNr,positionNr,:))' ...
            ,results.omega0,results.samplingFrequency,results.lags);
        r1WithPerturbationTheory(orientationNr,positionNr) = ...
            calculateR1WithSpectralDensity(spectralDensity1W0 ...
            ,spectralDensity2W0,results.dipolDipolConstant);
    end
end
r1WithPerturbationTheory
r1 = mean(results.r1WithPerturbationTheory(:,:,:,1),3)


corrFunc1W0 = results.allCorrelationFunction1W0Saver;
corrFunc2W0 = results.allCorrelationFunction2W0Saver;

disp('-- First caclulate R1 then average R1 --')
for atomCounter = 1:results.atomCounter
    for orientationNr = 1:length(orientationAngles)
        for positionNr = 1:length(positionAngles)
            
            [spectralDensity1W0, spectralDensity2W0] = ...
                calculateSpectralDensities(squeeze(corrFunc1W0(atomCounter,orientationNr ...
                ,positionNr,:))',squeeze(corrFunc2W0(atomCounter,orientationNr,positionNr,:))' ...
                ,results.omega0,results.samplingFrequency,results.lags);
            r1WithPerturbationTheory(atomCounter,orientationNr,positionNr) = ...
                calculateR1WithSpectralDensity(spectralDensity1W0 ...
                ,spectralDensity2W0,results.dipolDipolConstant);
        end
    end
end

squeeze(mean(r1WithPerturbationTheory,1))
r1 = mean(results.r1WithPerturbationTheory(:,:,:,1),3)



