clc
clear all %#ok<CLALL>

%% configure the System
configuration = readConfigurationFile('config.conf');
addpath(genpath(configuration.path2Library));
effectiveInteractionMyelinFraction = ...
    configuration.effectiveInteractionMyelinFraction;

optimizedParametersFromPaperCase1 = ...
    [0.15 0.5 2.51 1.2 14.9/effectiveInteractionMyelinFraction 16.3 1];
optimizedParametersFromPaperCase2 = ...
    [0.145 0.985 3.04 1.1 10/effectiveInteractionMyelinFraction 12.6 ...
    0.966];
optimizedParametersFromPaperCase3 = ...
    [0.14 0.59 2.57 0.86 14.7 16.0 0.99];
optimizedParametersFromPaperCase4 = ...
    [0.136 1.14 2.73 0.84 11.2 13.3 0.83];

%% load lipid data
lipidDataFieldsToLoad = configuration.lipidDataFieldsToLoad;
path2LipidData = configuration.path2LipidData;
lipidFieldNamesArray = getFieldNamesArray(lipidDataFieldsToLoad);
lipidDataFields = loadFieldsFromMatFile(path2LipidData ...
    ,lipidFieldNamesArray);

allLipidR1Rates = lipidDataFields.(lipidFieldNamesArray(1));
lipidAtomCounter = lipidDataFields.(lipidFieldNamesArray(4));
allLipidR1Rates = allLipidR1Rates(:,:,1:lipidAtomCounter);
lipidOrientations = rad2deg(lipidDataFields.(lipidFieldNamesArray(2)));
lipidPositions = rad2deg(lipidDataFields.(lipidFieldNamesArray(3)));

%% load water data
waterDataFieldsToLoad = configuration.waterDataFieldsToLoad;
path2WaterData = configuration.path2WaterData;
waterFieldNamesArray = getFieldNamesArray(waterDataFieldsToLoad);
waterDataFields = loadFieldsFromMatFile(path2WaterData ...
    ,waterFieldNamesArray);

waterAtomCounter = waterDataFields.(waterFieldNamesArray(4));
allWaterR1Rates = waterDataFields.(waterFieldNamesArray(1));
allWaterR1Rates = allWaterR1Rates(:,:,1:waterAtomCounter);
waterOrientations = rad2deg(waterDataFields.(waterFieldNamesArray(2)));
waterPositions = rad2deg(waterDataFields.(waterFieldNamesArray(3)));

lipidOrientationsCount = size(lipidOrientations,2);
waterOrientationsCount = size(waterOrientations,2);
if lipidOrientationsCount ~= waterOrientationsCount
    error(['Lipid Orientations (=' num2str(lipidOrientationsCount) ...
        ') and Water Orientations (=' num2str(waterOrientationsCount) ...
        ') are not the same.']);
else
    orientationsCount = lipidOrientationsCount;
end


%% Percentile
% TODO: calculate percentiles of the data set

%% calculate R1 rates

whichCase = configuration.whichCase;
switch whichCase
    case 1
        optimizedParameters = optimizedParametersFromPaperCase1;
        effectiveLipidR1Rates = mean(mean(allLipidR1Rates,3),2);
        effectiveWaterR1Rates = mean(mean(allWaterR1Rates,3),2);
        effectiveInteractionMyelinFraction = 0.32;
    case 2
        optimizedParameters = optimizedParametersFromPaperCase2;
        effectiveLipidR1Rates = median(median(allLipidR1Rates,3),2);
        effectiveWaterR1Rates = median(median(allWaterR1Rates,3),2);
        effectiveInteractionMyelinFraction = 0.32;
    case 3
        optimizedParameters = optimizedParametersFromPaperCase3;
        effectiveLipidR1Rates = mean(mean(allLipidR1Rates,3),2);
        effectiveWaterR1Rates = mean(mean(allWaterR1Rates,3),2);
        effectiveInteractionMyelinFraction = 1;
    case 4
        optimizedParameters = optimizedParametersFromPaperCase4;
        effectiveLipidR1Rates = median(median(allLipidR1Rates,3),2);
        effectiveWaterR1Rates = median(median(allWaterR1Rates,3),2);
        effectiveInteractionMyelinFraction = 1;
    otherwise 
        warning('Unknown Case.')
end

%% calculate the shifts

shapeFactor = optimizedParameters(7);
effectiveLipidR1RateShift = shapeFactor*( ...
    effectiveLipidR1Rates - effectiveLipidR1Rates(1));
effectiveWaterR1RateShift = shapeFactor*( ...
    effectiveWaterR1Rates - effectiveWaterR1Rates(1));

%% configure simulation
simulationTime = configuration.simulationTime;
timeSteps = configuration.timeSteps;
waterLipidHydrogenDensityRatio = configuration ...
    .waterLipidHydrogenDensityRatio;

timeAxis = linspace(0,simulationTime,timeSteps);
deltaT = timeAxis(2)-timeAxis(1);

%% calculate missing parameters
myelinWaterFraction = optimizedParameters(1);
solidMyelinFraction = (myelinWaterFraction ...
    /waterLipidHydrogenDensityRatio-myelinWaterFraction) ...
    *effectiveInteractionMyelinFraction;
freeWaterFraction = (1-myelinWaterFraction);

myelinWaterR1Offset = optimizedParameters(2);
solidMyelinR1Offset = optimizedParameters(3);

fittedFreeWaterR1Rates = optimizedParameters(4);
fittedMyelinWaterR1Rates = myelinWaterR1Offset+effectiveWaterR1RateShift;
fittedSolidMyelinR1RateRates = solidMyelinR1Offset+effectiveLipidR1RateShift;

exchangeRatesSM2MW = optimizedParameters(5);
exchangeRatesMW2FW = optimizedParameters(6);
exchangeRatesMW2SM = solidMyelinFraction/myelinWaterFraction ...
    *exchangeRatesSM2MW;
exchangeRatesFW2MW = myelinWaterFraction/freeWaterFraction ...
    *exchangeRatesMW2FW;

%% determination of R1 rate
predictedR1Rate = zeros(1,orientationsCount);

initialFreeWaterValue = -freeWaterFraction;
initialMyelinWaterValue = -myelinWaterFraction;
initialSolidMyelinValue = solidMyelinFraction;

for orientationNumber = 1:orientationsCount
    
    freeWaterPart = zeros(1,timeSteps);
    myelinWaterPart = zeros(1,timeSteps);
    solidMyelinPart = zeros(1,timeSteps);
    
    freeWaterPart(1) = initialFreeWaterValue;
    myelinWaterPart(1) = initialMyelinWaterValue;
    solidMyelinPart(1) = initialSolidMyelinValue;
    
    for timeStep = 1:timeSteps-1
        solidMyelinPart(timeStep+1) = solidMyelinPart(timeStep) ...
            +deltaT*(fittedSolidMyelinR1RateRates(orientationNumber) ...
            *(solidMyelinFraction-solidMyelinPart(timeStep)) ...
            -solidMyelinPart(timeStep)*exchangeRatesSM2MW ...
            +myelinWaterPart(timeStep)*exchangeRatesMW2SM);
        myelinWaterPart(timeStep+1) = myelinWaterPart(timeStep) ...
            +deltaT*(fittedMyelinWaterR1Rates(orientationNumber) ...
            *(myelinWaterFraction-myelinWaterPart(timeStep)) ...
            -exchangeRatesMW2SM*myelinWaterPart(timeStep) ...
            +exchangeRatesSM2MW*solidMyelinPart(timeStep) ...
            -exchangeRatesMW2FW*myelinWaterPart(timeStep) ...
            +exchangeRatesFW2MW*freeWaterPart(timeStep));
        freeWaterPart(timeStep+1) = freeWaterPart(timeStep) ...
            +deltaT*(fittedFreeWaterR1Rates ...
            *(freeWaterFraction-freeWaterPart(timeStep)) ...
            -exchangeRatesFW2MW*freeWaterPart(timeStep) ...
            +exchangeRatesMW2FW*myelinWaterPart(timeStep));
    end
end










 







