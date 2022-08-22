clc; clear all; close all;

%% system set up
saving = 1;
configuration = readConfigurationFile('orientationDependencyConfig.txt');
baseConfiguration = readConfigurationFile(configuration.baseConfig_path);
addpath(genpath(baseConfiguration.path2LibraryOnLocalMachine));
constants = readConstantsFile("../../constants.txt");

plotCorrelationFunctions = configuration.plotCorrelationFunctions;

dataToScale_paths.PLPC = configuration.corrFuncPLPCToScale_path;
dataToScale_paths.PSM = configuration.corrFuncPSMToScale_path;
dataToScale_paths.DOPS = configuration.corrFuncDOPSToScale_path;

scalingFactors_paths.PLPC = configuration.scalingFactorsPLPC_path;
scalingFactors_paths.PSM = configuration.scalingFactorsPSM_path;
scalingFactors_paths.DOPS = configuration.scalingFactorsDOPS_path;

fieldStrength = configuration.fieldStrength;
numberOfPhiAngles = configuration.numberOfPhiAngles;

lipidNames = fieldnames(dataToScale_paths);
for lipidNr = 1:length(lipidNames)
    lipidName = lipidNames{lipidNr};
    fprintf('Considered lipid: %s\n',lipidName);
    data_path = dataToScale_paths.(lipidName);
    scalFac_path = scalingFactors_paths.(lipidName);
    dataset = load(data_path);
    unscaledR1_path = calculateUnscaledR1ForOriDepPostProcessing(data_path, ...
        constants,fieldStrength);
    scaledR1_path = scaleUpR1WithScalingFactor(unscaledR1_path,scalFac_path);
    showResultsOfScaledUpR1(scaledR1_path);
    fprintf('\n\n');
    
end
