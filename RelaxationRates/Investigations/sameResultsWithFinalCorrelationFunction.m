clc
clear all
close all
% THAT does not work and I do not know why it is not working. Maybe I will
% check this later. There are one order higher predictions from this method
% compared to the results obtained from the simulation


plotting = false;

%% First calculate correlation function for whole system -> R1
configuration = readConfigurationFile('..\config.txt');
addpath(configuration.path2LibraryOnLocalMachine);

path2Data = "C:\Users\maxoe\Google Drive\Promotion\Results\performanceAnalysing\Server\20211008_Results_performanceAnalysing_Lipid_H_500ns_4ps_wh.mat";
%path2Data = "C:\Users\maxoe\Google Drive\Promotion\Results\performanceAnalysing\20211012_Results_performanceAnalysing_Lipid_H_500ns_1ps_nH40.mat";
data = load(path2Data);
dataConfiguration = data.configuration;
path2Conclusions = [configuration.path2ResultsOnLocalMachine ...
    dataConfiguration.kindOfResults '\Conclusions\'];
atomCounter = data.atomCounter;

correlationFunction0W0 = data.correlationFunction0W0Saver;
correlationFunction1W0 = data.correlationFunction1W0Saver;
correlationFunction2W0 = data.correlationFunction2W0Saver;

dt = data.deltaT;

orientationAngles = rad2deg(data.orientationAngles);
positionAngles = rad2deg(data.positionAngles);

r1Results = data.r1WithPerturbationTheory;
averagedR1Results = mean(r1Results(:,:,1:atomCounter),3);
dipolDipolConstant = data.dipolDipolConstant;
omega0 = data.constants.gyromagneticRatioOfHydrogenAtom ...
    *data.configuration.mainMagneticField;

for orientationNumber = 1:length(orientationAngles)
    for positionsNumber = 1:length(positionAngles)
        [spectralDensity1W0(orientationNumber,positionsNumber) ...
            ,spectralDensity2W0(orientationNumber,positionsNumber)] = ...
            calculateSpectralDensities( ...
            correlationFunction1W0(orientationNumber,positionsNumber) ...
            ,correlationFunction2W0(orientationNumber,positionsNumber) ...
            ,omega0,dt,data.lags); %#ok<SAGROW>
        r1WithCorrelationFunction(orientationNumber,positionsNumber) ...
            = calculateR1WithSpectralDensity( ...
            spectralDensity1W0(orientationNumber,positionsNumber) ...
            ,spectralDensity2W0(orientationNumber,positionsNumber) ...
            ,dipolDipolConstant); %#ok<SAGROW>
    end
end

averagedR1Results
r1WithCorrelationFunction