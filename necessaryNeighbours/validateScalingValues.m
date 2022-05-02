clc

configuration = readConfigurationFile('config.txt');
scalingRateConfig = 'scalingRateConfig.txt';
if configuration.runOnServer
    baseConfiguration = readConfigurationFile( ...
        configuration.path2BaseConfigurationOnServer);
    path2Results = baseConfiguration.path2ResultsOnServer;
    addpath(genpath(baseConfiguration.path2LibraryOnServer));
else
    baseConfiguration = readConfigurationFile( ...
        configuration.path2BaseConfigurationOnLocalMachine);
%     path2Results = baseConfiguration.path2ResultsOnLocalMachine;
    path2Results = "C:\Users\maxoe\Google Drive\Promotion\Results\LIPIDS\DOPS\nearestNeighbours\20220207_Results_relevantNearestNeighbours_DOPSlipid.mat";
    addpath(genpath(baseConfiguration.path2LibraryOnLocalMachine));
end

[~,~,~,path2LogFile] = ...
    setUpDependenciesBasedOnConfiguration(configuration);
deleteLogFile(path2LogFile);
logMessage('Starting Script to validate scaling rates.' ...
    ,path2LogFile);
printLineBreakToLogFile(path2LogFile);

%% Finding the right nearest neighbour cases
scalingRateConfiguration = readConfigurationFile(scalingRateConfig);
path2OldResultsFile = [path2Results ...
    scalingRateConfiguration.resultsFileName '.mat'];
path2OldResultsFile = path2Results;
results = load(path2OldResultsFile);
logMessage(sprintf('Determine scaling rates from results: %s'...
    ,path2OldResultsFile),path2LogFile,false);

[orientationsCount,positionsCount,calculatedAtoms,nNCasesCount] ...
    = size(results.r1WithPerturbationTheory);

averagedR1 = squeeze(mean(results.r1WithPerturbationTheory,3));
summedRates = zeros(1,nNCasesCount);


for referenceOrientationNr = 1:orientationsCount
    for referencePositionNr = 1:positionsCount        
        referenceR1 = squeeze(averagedR1(referenceOrientationNr ...
            ,referencePositionNr,:))';
        highestNNR1 = referenceR1(1);
        rates = highestNNR1./referenceR1;
        summedRates = summedRates + rates;
    end
end

averagedScalingRates = summedRates/(orientationsCount*positionsCount);
nearestNeighbourCasesOldSimulation = getValuesFromStringEnumeration( ...
    results.configuration.nearestNeighbourCases,';','numeric');

logMessage('Found the following cases with their scaling factors:' ...
    ,path2LogFile,false);
logMessage(sprintf('%.0f \t',nearestNeighbourCasesOldSimulation) ...
    ,path2LogFile,false);
logMessage(sprintf('%.6f\t',averagedScalingRates),path2LogFile,false);    


nearestNeighboursForValidation = ...
    getValuesFromStringEnumeration( ...
    scalingRateConfiguration.nearestNeighbourCases,';','numeric');
scalingFactors = zeros(1,length(nearestNeighboursForValidation));

for nNCaseNr = 1:length(nearestNeighboursForValidation)
    nNCase = nearestNeighboursForValidation(nNCaseNr);
    index = find(nearestNeighbourCasesOldSimulation == nNCase);
    scalingFactors(nNCaseNr) = averagedScalingRates(index);
end

logMessage(sprintf(['From which the following factors are used '...
    'according to: %s'],scalingRateConfig),path2LogFile,false);
logMessage(sprintf('%.0f \t',nearestNeighboursForValidation) ...
    ,path2LogFile,false);
logMessage(sprintf('%.3f\t',scalingFactors),path2LogFile,false);
printLineBreakToLogFile(path2LogFile);

clearvars -except scalingFactors nearestNeighboursForValidation ...
    resultsConfiguration configuration scalingRateConfiguration ...
    trajectoryX trajectoryY trajectoryZ simulationConfiguration ...
    averagedScalingRates nearestNeighbourCasesOldSimulation

[path2Data,path2Save,path2ConstantsFile,path2LogFile] = ...
    setUpDependenciesBasedOnConfiguration(configuration);

%% Analyse file name for simulation
fileName = configuration.fileName;
simulationDate = getSimulationDateFromFileName(fileName);
whichLipid = getLipidNameFromFileName(fileName);
waterModel = getWaterModelFromFileName(fileName);
formOfLayer = getFormOfLayerFromFileName(fileName);
waterMoleculesCount = getWaterMoleculesCountFromFileName(fileName);
constituent = getConstituentFromFileName(fileName);
composingMode = getComposingModeFromFileName(fileName);
samplingFrequency = getSamplingFrequencyFromFileName(fileName);
simTime = getSimulationTimeFromFileName(fileName);

path2File = sprintf('%s%s/%s.mat',path2Data,whichLipid,fileName);

logMessage(sprintf(['Data is simulated with the following ' ...
    'information: \n' ...
    '    GROMACS Simulation Date: %s \n' ...
    '    Lipid that is simulated: %s \n' ...
    '    Used Water Model: %s \n' ...
    '    Layer Form: %s \n' ...
    '    Number of Water Molecules: %s \n' ...
    '    Constituent of Simulation: %s \n' ...
    '    Composing mode in postprocessing: %s \n' ...
    '    Sampling Frequency: %s ps \n' ...
    '    Simulation Time: %s ns\n'],simulationDate,whichLipid,waterModel ...
    ,formOfLayer,waterMoleculesCount,constituent,composingMode ...
    ,samplingFrequency,simTime),path2LogFile,false);

logMessage(sprintf(['The following directories are used: \n' ...
    '    Path of Data: %s \n' ...
    '    Path of File: %s \n' ...
    '    Path where data is saved: %s \n' ...
    '    Path to Log File: %s \n'],path2Data,path2File,path2Save ...
    ,path2LogFile),path2LogFile,false);
logMessage('System was set up.',path2LogFile);

%% Load data
if ~configuration.dataLoaded
    logMessage('Start loading data.', path2LogFile);
    [trajectoryX,trajectoryY,trajectoryZ,simulationConfiguration] = ...
        loadTrajectoriesAndSimConfig(path2File);
    logMessage('Loading data finished',path2LogFile)
else
    logMessage('Data was already loaded in run before.',path2LogFile)
end

logMessage(sprintf('     The system has %i hydrogen atoms.' ...
    ,size(trajectoryX,1)),path2LogFile,false);

%% Define constants
logMessage('Defining constants.',path2LogFile,false);
constants = readConstantsFile(path2ConstantsFile);

dipolDipolConstant = 3/4*(constants.vaccumPermeability/(4*pi) ...
    *constants.hbar*constants.gyromagneticRatioOfHydrogenAtom^2)^2 ...
    /(constants.nanoMeter^6);
omega0 = constants.gyromagneticRatioOfHydrogenAtom ...
    *configuration.mainMagneticField;
samplingFrequency = samplingFrequency * constants.picoSecond;

%% Define simulation parameters
logMessage('Defining simulation parameters.',path2LogFile,false);

[numberOfHs,timeSteps] = size(trajectoryX);
logMessage(sprintf(['    Found %d hydrogen atoms at %d time steps of ' ...
    '%.3d s'],numberOfHs,timeSteps,samplingFrequency),path2LogFile,false);

lags = round(configuration.fractionForLags*timeSteps);
logMessage(sprintf(['    The lag is set to %d time steps, resulting ' ...
    'in a simulation time of %d s. Functionality comes into play ' ...
    'in correlation function'], lags, lags*samplingFrequency),path2LogFile ...
    ,false);

nearestNeighbourCases = nearestNeighboursForValidation;
nearestNeighbours = nearestNeighbourCases(1);
if nearestNeighbours >= numberOfHs
    logMessage(['The number of nearest neighbours is higher than the '...
        'number of possible atoms. PLEASE CHECK YOUR CONFIG FILE!'] ...
        ,path2LogFile);
    error(['The number of nearest neighbours is higher than the number' ...
        'of possible atoms. Please check your config file!']);
end
logMessage(sprintf(['    Analysing %.f nearst neighbours of overall ' ...
    '%.f hydrogen atoms'],nearestNeighbours,numberOfHs),path2LogFile ...
    ,false);

atomsToCalculate = configuration.atomsToCalculate;
startDateOfSimulation = datestr(now,'yyyymmdd');

orientationAngles = deg2rad(getValuesFromStringEnumeration( ...
    scalingRateConfiguration.fibreOrientations,';','numeric'));
fibreOrientationsCount = size(orientationAngles,2);
logMessage(['    Found the orientations' sprintf(' %.f',rad2deg( ...
    orientationAngles))],path2LogFile,false);

positionAngles = deg2rad(getValuesFromStringEnumeration( ...
    scalingRateConfiguration.myelinPositions,';','numeric'));
positionsInMyelinCount = size(positionAngles,2);
logMessage(['    Found the positions' sprintf(' %.f',rad2deg( ...
    positionAngles))],path2LogFile,false);

%% Preallocate some arrays
logMessage('Preallocation of some arrays.',path2LogFile,false);
nNCasesCount = length(nearestNeighboursForValidation);
r1WithPerturbationTheory = zeros(fibreOrientationsCount ...
    ,positionsInMyelinCount,atomsToCalculate,nNCasesCount);
scaledR1Rates = zeros(fibreOrientationsCount ...
    ,positionsInMyelinCount,atomsToCalculate,nNCasesCount);
meanPositions = single([mean(trajectoryX,2) mean(trajectoryY,2) ...
    mean(trajectoryZ,2)]);

atomCounter = 1;
atomIndex = zeros(1,atomsToCalculate);
calculationSteps = fibreOrientationsCount*positionsInMyelinCount;

randomSequenceOfAtoms = randperm(numberOfHs);

nearestNeighboursX = zeros(nearestNeighbours,timeSteps,'single');
nearestNeighboursY = zeros(nearestNeighbours,timeSteps,'single');
nearestNeighboursZ = zeros(nearestNeighbours,timeSteps,'single');
nearestNeighbourDistancesPow3 = zeros(nearestNeighbours,timeSteps ...
    ,'single');

polarAngle = zeros(nearestNeighbours,timeSteps,'single');
azimuthAngle = zeros(nearestNeighbours,timeSteps,'single');

sphericalHarmonicZerothOrder = zeros(nearestNeighbours,timeSteps ...
    ,'like',single(1j));
sphericalHarmonicFirstOrder = zeros(nearestNeighbours,timeSteps ...
    ,'like',single(1j));
sphericalHarmonicSecondOrder = zeros(nearestNeighbours,timeSteps ...
    ,'like',single(1j));

correlationFunctions0W0 = zeros(nNCasesCount,lags,'like',single(1j));
correlationFunctions1W0 = zeros(nNCasesCount,lags,'like',single(1j));
correlationFunctions2W0 = zeros(nNCasesCount,lags,'like',single(1j));

correlationFunction0W0Saver = zeros(fibreOrientationsCount ...
    ,positionsInMyelinCount,nNCasesCount,size(trajectoryX,2)...
    ,'like',single(1j));
correlationFunction1W0Saver = zeros(fibreOrientationsCount ...
    ,positionsInMyelinCount,nNCasesCount,size(trajectoryX,2) ...
    ,'like',single(1j));
correlationFunction2W0Saver = zeros(fibreOrientationsCount ...
    ,positionsInMyelinCount,nNCasesCount,size(trajectoryX,2) ...
    ,'like',single(1j));
logMessage('    Created correlation function saver.' ...
    , path2LogFile,false);

if configuration.reloadOldSimulation
    logMessage('Loading results from old simulation.',path2LogFile);
    oldSimulation = load(configuration.path2OldResults);
    checkForMismatchingConfigurations(configuration);
    logMessage(sprintf(['The old simulation was started on %s and ' ...
        'was last saved on %s.'],oldSimulation.startDateOfSimulation ...
        ,oldSimulation.lastSavingDate),path2LogFile,false);
    randomSequenceOfAtoms = oldSimulation.randomSequenceOfAtoms;
    atomCounter = oldSimulation.atomCounter;
    r1WithPerturbationTheory(:,:,1:atomCounter) = ...
        oldSimulation.r1WithPerturbationTheory(:,:,1:atomCounter);
    correlationFunction1W0Saver(:,:,1:atomCounter,:) = ...
        oldSimulation.correlationFunction1W0Saver(:,:,1:atomCounter,:);
    correlationFunction2W0Saver(:,:,1:atomCounter,:) = ...
        oldSimulation.correlationFunction2W0Saver(:,:,1:atomCounter,:);
end

printBreakLineToLogFile(path2LogFile);
atomTimer = [];

for atomNumber = randomSequenceOfAtoms(atomCounter:atomsToCalculate)
    atomTimerStart = tic;
    logMessage(sprintf('Selected atom number %i',atomNumber),path2LogFile);
    atomIndex(atomCounter) = atomNumber;
    
    inverseRotationMatrixPosition = eye(3);
    inverseRotationMatrixOrientation = eye(3);
    
    logMessage('    Relative positions.',path2LogFile,false);
    [trajectoryX,trajectoryY,trajectoryZ] ...
        = calculateRelativePositions(trajectoryX ...
        ,trajectoryY,trajectoryZ,atomNumber);
    
    [nearestNeighboursX,nearestNeighboursY,nearestNeighboursZ ...
        ,nearestNeighbourDistancesPow3] = findNearestNeighbours( ...
        nearestNeighbourCases(1),atomNumber,trajectoryX,trajectoryY ...
        ,trajectoryZ);
    
    for positionNumber = 1:positionsInMyelinCount
        positionTimer = tic;
        positionAngle = positionAngles(positionNumber);
        zAxis = [0 0 1];
        rotationMatrixPosition = get3DRotationMatrix( ...
            positionAngle,zAxis);
        for orientationNumber = 1:fibreOrientationsCount
            transformationTimer = tic;
            orientationAngle = orientationAngles(orientationNumber);
            logMessage(sprintf('=> Orientation: %i, position: %i' ...
                ,rad2deg(orientationAngle),rad2deg(positionAngle)) ...
                ,path2LogFile);
            yAxis = [0 1 0];
            rotationMatrixOrientation = get3DRotationMatrix( ...
                orientationAngle,yAxis);
            totalRotationMatrix = ...
                rotationMatrixOrientation*rotationMatrixPosition ...
                *inverseRotationMatrixPosition ...
                *inverseRotationMatrixOrientation ;
            logMessage('    Coordinate transformation.',path2LogFile ...
                ,false);
            [nearestNeighboursX,nearestNeighboursY,nearestNeighboursZ]  ...
                = rotateTrajectoriesWithRotationMatrix( ...
                totalRotationMatrix,nearestNeighboursX ...
                ,nearestNeighboursY,nearestNeighboursZ);
            
            logMessage('    Spherical coordinates.',path2LogFile,false);
            [polarAngle,azimuthAngle] = ...
                transformToSphericalCoordinates(nearestNeighboursX ...
                ,nearestNeighboursY,nearestNeighboursZ);
            
            logMessage('    Spherical harmonics.',path2LogFile,false);
            [sphericalHarmonicZerothOrder,sphericalHarmonicFirstOrder ...
                ,sphericalHarmonicSecondOrder] ...
                = calculateSphericalHarmonics(polarAngle,azimuthAngle ...
                ,nearestNeighbourDistancesPow3);
            polarAngle = zeros(nearestNeighbours,timeSteps,'single');
            azimuthAngle = zeros(nearestNeighbours,timeSteps,'single');
            
            logMessage('    Correlation function 0.',path2LogFile,false);
            correlationFunctions0W0 ...
                = calculateCorrelationFunctionsForMultipleNNCases( ...
                sphericalHarmonicZerothOrder,lags ...
                ,nearestNeighboursForValidation);
            correlationFunction0W0Saver(orientationNumber ...
                ,positionNumber,:,:) = (squeeze(correlationFunction0W0Saver( ...
                orientationNumber,positionNumber,:,:))*(atomCounter-1) + ...
                correlationFunctions0W0)/atomCounter;
            
            logMessage('    Correlation function 1.',path2LogFile,false);
            correlationFunctions1W0 ...
                = calculateCorrelationFunctionsForMultipleNNCases( ...
                sphericalHarmonicFirstOrder,lags ...
                ,nearestNeighboursForValidation);
            correlationFunction1W0Saver(orientationNumber ...
                ,positionNumber,:,:) = (squeeze(correlationFunction1W0Saver( ...
                orientationNumber,positionNumber,:,:))*(atomCounter-1) + ...
                correlationFunctions1W0)/atomCounter;
            
            logMessage('    Correlation function 2.',path2LogFile,false);
            correlationFunctions2W0 ...
                = calculateCorrelationFunctionsForMultipleNNCases( ...
                sphericalHarmonicSecondOrder,lags ...
                ,nearestNeighboursForValidation);
            correlationFunction2W0Saver(orientationNumber ...
                ,positionNumber,:,:) = (squeeze(correlationFunction2W0Saver( ...
                orientationNumber,positionNumber,:,:))*(atomCounter-1) + ...
                correlationFunctions2W0)/atomCounter;
            
            sphericalHarmonicZerothOrder = zeros(nearestNeighbours ...
                ,timeSteps,'like',single(1j));
            sphericalHarmonicFirstOrder = zeros(nearestNeighbours ...
                ,timeSteps,'like',single(1j));
            sphericalHarmonicSecondOrder = zeros(nearestNeighbours ...
                ,timeSteps,'like',single(1j));
            
            logMessage('    Spectral density and relaxation rates' ...
                ,path2LogFile,false);
            
            r1RatesMessage = "";
            scaledR1RatesMessage = "";
            for caseCounter = 1:nNCasesCount
                [spectralDensityW0,spectralDensity2W0] = ...
                    calculateSpectralDensities(correlationFunctions1W0( ...
                    caseCounter,:),correlationFunctions2W0( ...
                    caseCounter,:),omega0,samplingFrequency,lags);
                r1WithPerturbationTheory(orientationNumber,positionNumber ...
                    ,atomCounter,caseCounter) = calculateR1WithSpectralDensity( ...
                    spectralDensityW0,spectralDensity2W0,dipolDipolConstant);
                scaledR1Rates(orientationNumber,positionNumber ...
                    ,atomCounter,caseCounter) = ...
                    r1WithPerturbationTheory(orientationNumber ...
                    ,positionNumber,atomCounter,caseCounter) ...
                    *scalingFactors(caseCounter);
                r1RatesMessage = sprintf('%s  %i:%.5f',r1RatesMessage ...
                    ,nearestNeighboursForValidation(caseCounter) ...
                    ,r1WithPerturbationTheory(orientationNumber ...
                    ,positionNumber,atomCounter,caseCounter));
                scaledR1RatesMessage = sprintf('%s  %i:%.5f' ...
                    ,scaledR1RatesMessage ...
                    ,nearestNeighboursForValidation(caseCounter) ...
                    ,scaledR1Rates(orientationNumber,positionNumber ...
                    ,atomCounter,caseCounter));
            end

            logMessage(sprintf('Unscaled R1: %s',r1RatesMessage) ...
                ,path2LogFile);
            logMessage(sprintf('Scaled R1:   %s',scaledR1RatesMessage) ...
                ,path2LogFile);
            inverseRotationMatrixOrientation = inv( ...
                rotationMatrixOrientation);
        end
        inverseRotationMatrixPosition = inv(rotationMatrixPosition);
        logMessage(sprintf(['Finished position %i. Needed '...
            'time %.4f'],rad2deg(positionAngle),toc(positionTimer)) ...
            ,path2LogFile,false);
    end
    
    if mod(atomCounter,5) == 0
        lastSavingDate = datestr(now,'yyyymmdd_HHMM');
            createDataSavingObject();
            save(path2Save,'-struct','dataSavingObject','-v7.3');
            logMessage('Saved data',path2LogFile);
    end
    logMessage(sprintf('Calculated %i atoms',atomCounter),path2LogFile);
    atomTimer(end+1) = toc(atomTimerStart); %#ok<SAGROW>
    logMessage(sprintf([' ---> Average time for one atom: %s \n' ...
	'      Average time for one position/orientation: %s \n'...
	'      Approximately ready on: %s.'] ...
	,datestr(seconds(mean(atomTimer)),'HH:MM:SS') ...
	,datestr(seconds(mean(atomTimer)/calculationSteps),'HH:MM:SS') ...
	,datetime('now')+seconds(mean(atomTimer) ...
	*(atomsToCalculate-atomCounter))),path2LogFile,false);
    printDottedBreakLineToLogFile(path2LogFile);
    atomCounter = atomCounter + 1;
end


