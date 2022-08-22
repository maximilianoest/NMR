clc; clear all; close all;

%% system set up
saving = 1;
nNConfiguration = readConfigurationFile('nNConfig.txt');
baseConfiguration =  readConfigurationFile(nNConfiguration.baseConfig_path);
addpath(genpath(baseConfiguration.path2LibraryOnLocalMachine));

plotCorrelationFunctions = nNConfiguration.plotCorrelationFunctions;

offsetSuppressionRegion = getValuesFromStringEnumeration( ...
    nNConfiguration.offsetSuppressionRegion,';','numeric');
cutOffFraction = nNConfiguration.cutOffFraction;
fieldStrengths = getValuesFromStringEnumeration( ...
    nNConfiguration.fieldStrengths,';','numeric');

corrFuncForScalFac_paths.PLPC = nNConfiguration.corrFuncPLPCForScalFac_path;
corrFuncForScalFac_paths.PSM = nNConfiguration.corrFuncPSMForScalFac_path;
corrFuncForScalFac_paths.DOPS = nNConfiguration.corrFuncDOPSForScalFac_path;

corrFuncToScale_paths.PLPC = nNConfiguration.corrFuncPLPCToScale_path;
corrFuncToScale_paths.PSM = nNConfiguration.corrFuncPSMToScale_path;
corrFuncToScale_paths.DOPS = nNConfiguration.corrFuncDOPSToScale_path;


lipidNames = fieldnames(corrFuncForScalFac_paths);
for fieldStrengthNr = 1:size(fieldStrengths,2)
    fieldStrength = fieldStrengths(fieldStrengthNr);
    for lipidNr = 1:size(lipidNames,1)
        lipidName = lipidNames{lipidNr};
        corrFunc_path = corrFuncForScalFac_paths.(lipidName);
        corrFuncValidation_path = ...
            corrFuncToScale_paths.(lipidName);
        if ~exist(corrFunc_path,'file')
            warning('nNPostProcessing:fileForLipidNotFound', ...
                'The File: \n %s \ncannot be found. Skip to next file.', ...
                corrFunc_path);
            continue
        end
        %% determine R1 from corr func
        if fieldStrengthNr == 1
            if plotCorrelationFunctions ~= 0
                plotCorrelationFunctions = 1;
            end
        else
            plotCorrelationFunctions = 0;
        end
        r1ForDifferentNNForScalFac_path = calculateR1ForNNPostProcessing( ...
            corrFunc_path,fieldStrength,offsetSuppressionRegion, ...
            cutOffFraction,plotCorrelationFunctions,'');
        
        %% determine scaling factors
        scalFac_path = determineScalingFactors( ...
            r1ForDifferentNNForScalFac_path,fieldStrength);
        
        %% scaling up R1 with less NN from other data set
        if fieldStrengthNr == 1
            if plotCorrelationFunctions ~= 0
                plotCorrelationFunctions = 1;
            end
        else
            plotCorrelationFunctions = 0;
        end
        r1ForDifferentNNToScale_path = ...
            calculateR1ForNNPostProcessing(corrFuncValidation_path, ...
            fieldStrength,offsetSuppressionRegion,cutOffFraction, ...
            plotCorrelationFunctions,'Validation');
        
        r1ScaledUp_path = scaleUpR1WithScalingFactorsForValidation( ...
            r1ForDifferentNNToScale_path,scalFac_path);
        
        showThetaDependenceOfR1ForTwoDataSets(...
            r1ForDifferentNNForScalFac_path, ...
            r1ForDifferentNNToScale_path,saving);
        %%
        
        compareDistributionsOfDatasets(corrFuncValidation_path,saving);
        % show all distributions in one plot
        % (corrFuncValidation_path,corrFunc_path,saving)
        showDifferentDistributionsInOnePLot(corrFuncValidation_path, ...
            corrFunc_path,'X',0.15,1)
        
        %% compare to R1 calculated with higher NN
        compareNNDependentScaledR1(r1ScaledUp_path,saving);
        
        
        
    end
    
end


