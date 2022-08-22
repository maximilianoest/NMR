 %% Set Up System

clc
clearvars -except trajectoryX trajectoryY trajectoryZ simulationConfiguration
configuration = readConfigurationFile('config.txt');
if configuration.runOnServer
    baseConfiguration = readConfigurationFile( ...
        configuration.path2BaseConfigurationOnServer);
    addpath(genpath(baseConfiguration.path2LibraryOnServer));
else
    baseConfiguration = readConfigurationFile( ...
        configuration.path2BaseConfigurationOnLocalMachine);
    addpath(genpath(baseConfiguration.path2LibraryOnLocalMachine))
end

[path2Data,path2Save,path2ConstantsFile,path2LogFile] = ...
    setUpDependenciesBasedOnConfiguration(configuration);
deleteLogFile(path2LogFile);
logMessage('Starting Script to analyse nearest Neighbours dependence.' ...
    ,path2LogFile);

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
if ~exist(path2File,'file')
    logMessage(sprintf('FILE: %s CANNOT BE FOUND!',path2File) ...
        ,path2LogFile, false);
    error('File cannot be found. See log file.');
else
    logMessage(sprintf('File: %s exists and will be loaded.',path2File) ...
        ,path2LogFile, false);
end

logMessage(sprintf(['Data are simulated with the following ' ...
    'information: \n' ...
    '    GROMACS Simulation Date: %s \n' ...
    '    Lipid that is simulated: %s \n' ...
    '    Used Water Model: %s \n' ...
    '    Layer Form: %s \n' ...
    '    Number of Water Molecules: %s \n' ...
    '    Constituent of Simulation: %s \n' ...
    '    Composing mode in postprocessing: %s \n' ...
    '    Sampling Frequency: %f ps \n' ...
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

nearestNeighbourCases = getValuesFromStringEnumeration( ...
    configuration.nearestNeighbourCases,';','numeric');
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
    configuration.fibreOrientations,';','numeric'));
fibreOrientationsCount = size(orientationAngles,2);
logMessage(['    Found the orientations' sprintf(' %.f',rad2deg( ...
    orientationAngles))],path2LogFile,false);

positionAngles = deg2rad(getValuesFromStringEnumeration( ...
    configuration.myelinPositions,';','numeric'));
positionsInMyelinCount = size(positionAngles,2);
logMessage(['    Found the positions' sprintf(' %.f',rad2deg( ...
    positionAngles))],path2LogFile,false);

correlationFunction0W0Saver = zeros(length(nearestNeighbourCases) ...
    ,fibreOrientationsCount,positionsInMyelinCount,size(trajectoryX,2) ...
    ,'like',single(1j));
correlationFunction1W0Saver = zeros(length(nearestNeighbourCases) ...
    ,fibreOrientationsCount,positionsInMyelinCount,size(trajectoryX,2) ...
    ,'like',single(1j));
correlationFunction2W0Saver = zeros(length(nearestNeighbourCases) ...
    ,fibreOrientationsCount,positionsInMyelinCount,size(trajectoryX,2) ...
    ,'like',single(1j));

sumCorrelationFunction0W0Saver = zeros(length(nearestNeighbourCases) ...
    ,fibreOrientationsCount,positionsInMyelinCount,size(trajectoryX,2) ...
    ,'like',single(1j));
sumCorrelationFunction1W0Saver = zeros(length(nearestNeighbourCases) ...
    ,fibreOrientationsCount,positionsInMyelinCount,size(trajectoryX,2) ...
    ,'like',single(1j));
sumCorrelationFunction2W0Saver = zeros(length(nearestNeighbourCases) ...
    ,fibreOrientationsCount,positionsInMyelinCount,size(trajectoryX,2) ...
    ,'like',single(1j));

logMessage('    Created correlation function saver.' ...
    , path2LogFile,false);

%% Start simulation
logMessage('Preallocation of some other arrays.',path2LogFile,false);

r1WithPerturbationTheory = zeros(fibreOrientationsCount ...
    ,positionsInMyelinCount,atomsToCalculate);
meanPositions = single([mean(trajectoryX,2) mean(trajectoryY,2) ...
    mean(trajectoryZ,2)]);
inverseRotationMatrixPosition = eye(3);
inverseRotationMatrixOrientation = eye(3);

atomCounter = 1;
atomIndex = zeros(1,atomsToCalculate);
calculationSteps = fibreOrientationsCount*positionsInMyelinCount;

randomSequenceOfAtoms = reloadOldAtomSequenceToAvoidOverlappingDataSets();

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

correlationFunction0W0 = zeros(nearestNeighbourCases(1),lags,'like' ...
    ,single(1j));
correlationFunction1W0 = zeros(nearestNeighbourCases(1),lags,'like' ...
    ,single(1j));
correlationFunction2W0 = zeros(nearestNeighbourCases(1),lags,'like' ...
    ,single(1j));

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
    sumCorrelationFunction1W0Saver(:,:,1:atomCounter,:) = ...
        oldSimulation.correlationFunction1W0Saver(:,:,1:atomCounter,:);
    sumCorrelationFunction2W0Saver(:,:,1:atomCounter,:) = ...
        oldSimulation.correlationFunction2W0Saver(:,:,1:atomCounter,:);
end

printBreakLineToLogFile(path2LogFile);
atomTimer = [];

for atomNumber = randomSequenceOfAtoms(atomCounter:atomsToCalculate)
    atomTimerStart = tic;
    logMessage(sprintf('Selected atom number %i',atomNumber),path2LogFile);
    atomIndex(atomCounter) = atomNumber;
    
    [trajectoryX,trajectoryY,trajectoryZ] ...
        = calculateRelativePositions(trajectoryX ...
        ,trajectoryY,trajectoryZ,atomNumber);
    
    [nearestNeighboursX,nearestNeighboursY,nearestNeighboursZ ...
        ,nearestNeighbourDistancesPow3] = findNearestNeighbours( ...
        nearestNeighbourCases(1),atomNumber,trajectoryX,trajectoryY ...
        ,trajectoryZ);
    
    for positionNumber = 1:positionsInMyelinCount
        positionAngle = positionAngles(positionNumber);
        zAxis = [0 0 1];
        rotationMatrixPosition = get3DRotationMatrix( ...
            positionAngle,zAxis);
        for orientationNumber = 1:fibreOrientationsCount
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

            [nearestNeighboursX,nearestNeighboursY,nearestNeighboursZ]  ...
                = rotateTrajectoriesWithRotationMatrix( ...
                totalRotationMatrix,nearestNeighboursX ...
                ,nearestNeighboursY,nearestNeighboursZ);
            
            [polarAngle,azimuthAngle] = ...
                transformToSphericalCoordinates(nearestNeighboursX ...
                ,nearestNeighboursY,nearestNeighboursZ);
            
            [sphericalHarmonicZerothOrder,sphericalHarmonicFirstOrder ...
                ,sphericalHarmonicSecondOrder] ...
                = calculateSphericalHarmonics(polarAngle,azimuthAngle ...
                ,nearestNeighbourDistancesPow3);
            polarAngle = zeros(nearestNeighbours,timeSteps,'single');
            azimuthAngle = zeros(nearestNeighbours,timeSteps,'single');
            
            correlationFunction0W0 ...
                = calculateCorrelationFunctionForNearestNeighbours( ...
                sphericalHarmonicZerothOrder,lags);
            correlationFunction1W0 ...
                = calculateCorrelationFunctionForNearestNeighbours( ...
                sphericalHarmonicFirstOrder,lags);
            correlationFunction2W0 ...
                = calculateCorrelationFunctionForNearestNeighbours( ...
                sphericalHarmonicSecondOrder,lags);
            
            sphericalHarmonicZerothOrder = zeros(nearestNeighbours ...
                ,timeSteps,'like',single(1j));
            sphericalHarmonicFirstOrder = zeros(nearestNeighbours ...
                ,timeSteps,'like',single(1j));
            sphericalHarmonicSecondOrder = zeros(nearestNeighbours ...
                ,timeSteps,'like',single(1j));
            
            r1Rates = "";
            for caseCounter = 1:length(nearestNeighbourCases)
                nearestNeighboursReduction = ...
                    nearestNeighbourCases(caseCounter);
                reducedCorrelationFunction0W0 = sum( ...
                    correlationFunction0W0( ...
                    1:nearestNeighboursReduction,:),1);
                reducedCorrelationFunction1W0 = sum( ...
                    correlationFunction1W0( ...
                    1:nearestNeighboursReduction,:),1);
                reducedCorrelationFunction2W0 = sum( ...
                    correlationFunction2W0( ...
                    1:nearestNeighboursReduction,:),1);
                
                correlationFunction0W0Saver(caseCounter,orientationNumber ...
                    ,positionNumber,:) = (squeeze( ...
                    correlationFunction0W0Saver(caseCounter,orientationNumber ...
                    ,positionNumber,:))'*(atomCounter-1) + ...
                    reducedCorrelationFunction0W0)/atomCounter;
                correlationFunction1W0Saver(caseCounter,orientationNumber ...
                    ,positionNumber,:) = (squeeze( ...
                    correlationFunction1W0Saver(caseCounter,orientationNumber ...
                    ,positionNumber,:))'*(atomCounter-1) + ...
                    reducedCorrelationFunction1W0)/atomCounter;
                correlationFunction2W0Saver(caseCounter,orientationNumber ...
                    ,positionNumber,:) = (squeeze( ...
                    correlationFunction2W0Saver(caseCounter,orientationNumber ...
                    ,positionNumber,:))'*(atomCounter-1) + ...
                    reducedCorrelationFunction2W0)/atomCounter;
                
                sumCorrelationFunction0W0Saver(caseCounter,orientationNumber ...
                    ,positionNumber,:) = squeeze(sumCorrelationFunction0W0Saver( ...
                    caseCounter,orientationNumber,positionNumber,:))' ...
                    + reducedCorrelationFunction0W0;
                sumCorrelationFunction1W0Saver(caseCounter,orientationNumber ...
                    ,positionNumber,:) = squeeze(sumCorrelationFunction1W0Saver( ...
                    caseCounter,orientationNumber,positionNumber,:))' ...
                    + reducedCorrelationFunction1W0;
                sumCorrelationFunction2W0Saver(caseCounter,orientationNumber ...
                    ,positionNumber,:) = squeeze(sumCorrelationFunction2W0Saver( ...
                    caseCounter,orientationNumber,positionNumber,:))' ...
                    + reducedCorrelationFunction2W0;
                
                [spectralDensityW0,spectralDensity2W0] = ...
                    calculateSpectralDensities( ...
                    reducedCorrelationFunction1W0 ...
                    ,reducedCorrelationFunction2W0,omega0 ...
                    ,samplingFrequency,lags);
                
                r1WithPerturbationTheory(orientationNumber ...
                    ,positionNumber,atomCounter,caseCounter) = ...
                    calculateR1WithSpectralDensity(spectralDensityW0 ...
                    ,spectralDensity2W0,dipolDipolConstant);
                r1Rates = sprintf('%s  %i:%.5f',r1Rates ...
                    ,nearestNeighboursReduction ...
                    ,r1WithPerturbationTheory(orientationNumber ...
                    ,positionNumber,atomCounter,caseCounter));
            end 
            logMessage(r1Rates,path2LogFile);
            inverseRotationMatrixOrientation = inv( ...
                rotationMatrixOrientation);
        end
        inverseRotationMatrixPosition = inv(rotationMatrixPosition);
    end
    
    if mod(atomCounter,1) == 0
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

