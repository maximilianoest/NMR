clc
clear all %#ok<CLALL>

%% configure the System
configuration = readConfigurationFile('config.conf');
addpath(genpath(configuration.path2Library));

%% get data
lipidData = load([configuration.path2LipidData ...
    configuration.lipidFileName]);
waterData = load([configuration.path2WaterData ...
    configuration.waterFileName]);

%% load lipid data
allLipidR1Rates = lipidData.calculatedR1Rates;
lipidOrientations = rad2deg(lipidData.positionAngles);

%% load water data
allWaterR1Rates = waterData.calculatedR1Rates;
waterOrientations = rad2deg(waterData.orientationAngles);

%% check for same data
lipidOrientationsCount = size(lipidOrientations,2);
waterOrientationsCount = size(waterOrientations,2);
if lipidOrientationsCount == waterOrientationsCount
    orientationsCount = lipidOrientationsCount;
    orientations = lipidOrientations;
else
    error(['Lipid Orientations (=' num2str(lipidOrientationsCount) ...
        ') and Water Orientations (=' num2str(waterOrientationsCount) ...
        ') are not the same!']);
end

%% decide which data to calculate and show
if lipidData.B0 == waterData.B0
    whichCaseForOptimizedParameters = ...
        configuration.whichCaseForOptimizedParameters;
    optimizedParameters = getHardcodedOptimizedParametersFromPaper( ...
        whichCaseForOptimizedParameters,lipidData.B0);
    
    whichCaseForEffectiveR1Rates = ...
        configuration.whichCaseForEffectiveR1Rates;
    [effectiveLipidR1Rates,effectiveWaterR1Rates] = ...
        getEffectiveR1RatesForCase(whichCaseForEffectiveR1Rates ...
        ,lipidData,waterData);
else
    error('The magnetic fields do not have the same strengths!')
end


%% calculate the shifts
shapeFactor = optimizedParameters.shapeFactor;
effectiveLipidR1RateShift = shapeFactor*( ...
    effectiveLipidR1Rates - effectiveLipidR1Rates(1));
effectiveWaterR1RateShift = shapeFactor*( ...
    effectiveWaterR1Rates - effectiveWaterR1Rates(1));

%% configure simulation
simulationTime = configuration.simulationTime;
timeSteps = configuration.timeSteps;

timeAxis = linspace(0,simulationTime,timeSteps);
deltaT = timeAxis(2)-timeAxis(1);

%% calculate missing parameters
effectiveInteractionMyelinFraction = ...
    configuration.effectiveInteractionMyelinFraction;
waterLipidHydrogenDensityRatio = configuration ...
    .waterLipidHydrogenDensityRatio;

myelinWaterFraction = optimizedParameters.myelinWaterFraction;
solidMyelinFraction = (myelinWaterFraction ...
    /waterLipidHydrogenDensityRatio-myelinWaterFraction) ...
    *effectiveInteractionMyelinFraction;
freeWaterFraction = (1-myelinWaterFraction);

myelinWaterR1Offset = optimizedParameters.myelinWaterR1Offset;
fittedMyelinWaterR1Rates = myelinWaterR1Offset ...
    +effectiveWaterR1RateShift;

solidMyelinR1Offset = optimizedParameters.solidMyelinR1Offset;
fittedSolidMyelinR1RateRates = solidMyelinR1Offset ...
    +effectiveLipidR1RateShift;

freeWaterR1Offset = optimizedParameters.freeWaterR1Offset;

exchangeRatesSM2MW = optimizedParameters.exchangeRatesSM2MW;
exchangeRatesMW2FW = optimizedParameters.exchangeRatesMW2FW;
exchangeRatesMW2SM = solidMyelinFraction/myelinWaterFraction ...
    *exchangeRatesSM2MW;
exchangeRatesFW2MW = myelinWaterFraction/freeWaterFraction ...
    *exchangeRatesMW2FW;

%% calculate R1 rate
predictedR1Rates = zeros(1,orientationsCount);

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
            +deltaT*(freeWaterR1Offset ...
            *(freeWaterFraction-freeWaterPart(timeStep)) ...
            -exchangeRatesFW2MW*freeWaterPart(timeStep) ...
            +exchangeRatesMW2FW*myelinWaterPart(timeStep));
    end
    longitudinalMagnetizsation = myelinWaterPart+freeWaterPart;
    expontentialCurve = @(r1Fitted,data) (1-2*exp(-r1Fitted(1)*data));
    startValue = [ 1.2 ]; %#ok<NBRAK>
    options = optimset('Display','off');
    predictedR1Rates(orientationNumber) = ...
        lsqcurvefit(expontentialCurve,startValue,timeAxis ...
        ,longitudinalMagnetizsation,[],[],options);
end

%% plotting
informationText = {'Data dates:' ...
    ,['Lipid: ' lipidData.startDateOfSimulation] ...
    ,['Water: ' waterData.startDateOfSimulation]};
figs(1) = figure(1);
plot(orientations,predictedR1Rates,'LineWidth',1.5)
title([whichCaseForEffectiveR1Rates ' and ' ...
    whichCaseForOptimizedParameters])
xlabel('Angle \theta [°]')
ylabel('Relaxation Rate R_1 [Hz]')
grid on
annotation('textbox',[.9 .5 .1 .2],'String',informationText ...
    ,'EdgeColor','none')
saveas(figs,[configuration.path2SaveFigs whichCaseForEffectiveR1Rates ...
    '_' configuration.fileNameToSaveFigs])



 







