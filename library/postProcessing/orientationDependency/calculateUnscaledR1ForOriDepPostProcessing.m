function fullFilePath = calculateUnscaledR1ForOriDepPostProcessing( ...
    dataset_path,constants,fieldStrength)
% orientation dependency configuration
configuration = readConfigurationFile('orientationDependencyConfig.txt');

% constants stuff
dipolDipolConstant = 3/4*(constants.vaccumPermeability/(4*pi) ...
    *constants.hbar*constants.gyromagneticRatioOfHydrogenAtom^2)^2 ...
    /(constants.nanoMeter^6);
gyromagneticRatio = constants.gyromagneticRatioOfHydrogenAtom;

% contribution of magnetic field strength
omega0 = gyromagneticRatio*fieldStrength;

dataset = load(dataset_path);
analyseAndShowSimulationParametersOfResultsDataset(dataset);
whichLipid = dataset.whichLipid;
matlabSimulationDate = dataset.startDateOfSimulation;
nearestNeighbours = dataset.nearestNeighbours;
atomCounter = dataset.atomCounter;
gromacsSimulationDate = dataset.simulationDate;
gromacsFileName = dataset.fileName;

deltaT = dataset.samplingFrequency;
corrFunc1W0 = double(dataset.correlationFunction1W0Saver);
corrFunc2W0 = double(dataset.correlationFunction2W0Saver);
corrFuncLength = size(corrFunc1W0,3);

orientationAngles = rad2deg(dataset.orientationAngles);
positionAngles = rad2deg(dataset.positionAngles);

offSetSuppression = getValuesFromStringEnumeration( ...
    configuration.offsetSuppressionRegion,';','numeric');
offsetSuppressionRegion = round(offSetSuppression(1)*corrFuncLength) ...
    :round(offSetSuppression(2)*corrFuncLength);

remainingRegion = 1:round(corrFuncLength*configuration.cutOffFraction)-1;

if isfield(dataset.configuration, 'AAXX_WARNING')
    warning('Phi was changed to make the files coincident!')
end

for orientationNr = 1:length(orientationAngles)
    for positionNr = 1:length(positionAngles)
        
        % offset suppression due to static part of Hamiltonian
        offsetSuppressedCorrFuncs1W0(orientationNr,positionNr,:) = ...
            suppressOffsetOfCorrFunc(corrFunc1W0( ...
            orientationNr,positionNr,:), ...
            offsetSuppressionRegion); %#ok<AGROW>
        offsetSuppressedCorrFuncs2W0(orientationNr,positionNr,:) = ...
            suppressOffsetOfCorrFunc(corrFunc2W0( ...
            orientationNr,positionNr,:), ...
            offsetSuppressionRegion); %#ok<AGROW>
        
        % cut off part of correlation function to avaoid rise of
        % correlation function with increasing tau caused by FT to
        % determine correlation function.
        cuttedCorrFuncs1W0(orientationNr, ...
            positionNr,:) = offsetSuppressedCorrFuncs1W0( ...
            orientationNr,positionNr,remainingRegion); %#ok<AGROW>
        cuttedCorrFuncs2W0(orientationNr, ...
            positionNr,:) = offsetSuppressedCorrFuncs2W0( ...
            orientationNr,positionNr,remainingRegion); %#ok<AGROW>
    end
end

% calculate R1 from cutted and offset suppressed correlation function
lags = size(cuttedCorrFuncs1W0,3);
for orientationNr = 1:length(orientationAngles)
    for positionNr = 1:length(positionAngles)
        unscaledR1_theta_phi(orientationNr,positionNr) = ...
            calculateR1FromCorrFuncWithExplIntegration( ...
            squeeze(cuttedCorrFuncs1W0(orientationNr, ...
            positionNr,:))',squeeze(cuttedCorrFuncs2W0( ...
            orientationNr,positionNr,:))',deltaT,omega0,lags, ...
            dipolDipolConstant); %#ok<AGROW>
    end
end

savingPath = initializeSystemForSavingR1();
fieldstrengthString = strrep(num2str(fieldStrength),'.','');
savingName = sprintf('%s_%sTesla_%iH_%iNN_r1%s_unscaledRelaxationRatesR1', ...
    whichLipid,fieldstrengthString,atomCounter,nearestNeighbours, ...
    matlabSimulationDate);
fullFilePath = [savingPath savingName '.mat'];
save(fullFilePath,'unscaledR1_theta_phi','whichLipid','fieldStrength', ...
    'matlabSimulationDate','orientationAngles','positionAngles', ...
    'nearestNeighbours','atomCounter','gromacsSimulationDate', ...
    'gromacsFileName');

end
