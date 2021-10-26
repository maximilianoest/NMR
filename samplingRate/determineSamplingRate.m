%% Set Up System
clc
clear all %#ok<CLALL>

configuration = readConfigurationFile('config.txt');
if configuration.runOnServer
    addpath(genpath(configuration.path2LibraryOnServer));
else
    addpath(genpath(configuration.path2LibraryOnLocalMachine))
end

[path2Data,path2Save,path2ConstantsFile,path2LogFile] = ...
    setUpSystemBasedOnMachine(configuration);
deleteLogFile(path2LogFile);
fileName = configuration.fileName;

logMessage('System is set up.', path2LogFile);

%% Load data
logMessage('Start loading data.', path2LogFile);
[trajectoryX,trajectoryY,trajectoryZ] = loadTrajectoriesFromData( ...
    configuration,path2Data);
logMessage('Data successfully loaded.', path2LogFile);

%% Define constants
logMessage('Defining constants.',path2LogFile,false);
deltaT = configuration.deltaT;
constants = readConstantsFile(path2ConstantsFile);

dipolDipolConstant = 3/4*(constants.vaccumPermeability/(4*pi) ...
    *constants.hbar*constants.gyromagneticRatioOfHydrogenAtom^2)^2 ...
    /(constants.nanoMeter^6);
omega0 = constants.gyromagneticRatioOfHydrogenAtom ...
    *configuration.mainMagneticField; 

%% Define simulation parameters
logMessage('Defining simulation parameters.',path2LogFile,false);

[numberOfHs,timeSteps] = size(trajectoryX);
logMessage(sprintf(['    Found %d hydrogen atoms at %d time steps of ' ...
    '%.3d s'],numberOfHs,timeSteps,deltaT),path2LogFile,false);

lags = round(configuration.fractionForLags*timeSteps);
logMessage(sprintf(['    The lag is set to %d time steps, resulting ' ...
    'in a simulation time of %d s. NO FUNCTIONALITY'], lags ...
    ,lags*deltaT),path2LogFile,false);

nearestNeighbours = configuration.nearestNeighbours;
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

correlationFunction0W0Saver = zeros(fibreOrientationsCount ...
    ,positionsInMyelinCount,size(trajectoryX,2),'like',single(1j));
correlationFunction1W0Saver = zeros(fibreOrientationsCount ...
    ,positionsInMyelinCount,size(trajectoryX,2),'like',single(1j));
correlationFunction2W0Saver = zeros(fibreOrientationsCount ...
    ,positionsInMyelinCount,size(trajectoryX,2),'like',single(1j));
logMessage('    Created correlation function saver.' ...
    , path2LogFile,false);

%% Start simulation
logMessage('Preallocation of some other arrays.',path2LogFile,false);

r1WithPerturbationTheory = zeros(fibreOrientationsCount ...
    ,positionsInMyelinCount,atomsToCalculate);
meanPositions = single([mean(trajectoryX,2) mean(trajectoryY,2) ...
    mean(trajectoryZ,2)]);

atomCounter = 1;
atomIndex = zeros(1,atomsToCalculate);
calculationSteps = fibreOrientationsCount*positionsInMyelinCount;

randomSequenceOfAtoms = randperm(numberOfHs);

relativeX = zeros(numberOfHs,timeSteps,'like',single(1));
relativeY= zeros(numberOfHs,timeSteps,'like',single(1));
relativeZ = zeros(numberOfHs,timeSteps,'like',single(1));

nearestNeighboursX = zeros(nearestNeighbours,timeSteps,'like',single(1));
nearestNeighboursY = zeros(nearestNeighbours,timeSteps,'like',single(1));
nearestNeighboursZ = zeros(nearestNeighbours,timeSteps,'like',single(1));
nearestNeighbourDistancesPow3 = zeros(nearestNeighbours,timeSteps ...
    ,'like',single(1));

rotatedX = zeros(nearestNeighbours,timeSteps,'like',single(1));
rotatedY = zeros(nearestNeighbours,timeSteps,'like',single(1));
rotatedZ = zeros(nearestNeighbours,timeSteps,'like',single(1));

polarAngle = zeros(nearestNeighbours,timeSteps,'like',single(1));
azimuthAngle = zeros(nearestNeighbours,timeSteps,'like',single(1));

sphericalHarmonicZerothOrder = zeros(nearestNeighbours,timeSteps ...
    ,'like',single(1j));
sphericalHarmonicFirstOrder = zeros(nearestNeighbours,timeSteps ...
    ,'like',single(1j));
sphericalHarmonicSecondOrder = zeros(nearestNeighbours,timeSteps ...
    ,'like',single(1j));

correlationFunction0W0 = zeros(1,lags,'like',single(1j));
correlationFunction1W0 = zeros(1,lags,'like',single(1j));
correlationFunction2W0 = zeros(1,lags,'like',single(1j));

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
    timeTracks = oldSimulation.timeTracks;
end

logMessage('Set up done. Starting Simulation.',path2LogFile);
printBreakLineToLogFile(path2LogFile);

% for atomNumber = 1:numberOfHs
for atomNumber = randomSequenceOfAtoms(atomCounter:end)
    overallForAtom = tic;
    logMessage(sprintf('Selected atom number %i',atomNumber),path2LogFile);
    atomIndex(atomCounter) = atomNumber;
    
    logMessage('Calcualting relative positions.',path2LogFile,false);
    [relativeX,relativeY,relativeZ] ...
        = calculateRelativePositions(trajectoryX ...
        ,trajectoryY,trajectoryZ,atomNumber);
    
    logMessage('Nearest neighbours.',path2LogFile,false);
    [nearestNeighboursX,nearestNeighboursY,nearestNeighboursZ ...
        ,nearestNeighbourDistancesPow3] = findNearestNeighbours( ...
        nearestNeighbours,atomNumber,relativeX,relativeY ...
        ,relativeZ);
    
    for positionNumber = 1:positionsInMyelinCount
        positionAngle = positionAngles(positionNumber);
        zAxis = [0 0 1];
        rotationMatrixPosition = get3DRotationMatrix( ...
            positionAngle,zAxis);
        for orientationNumber = 1:fibreOrientationsCount
            orientationAngle = orientationAngles(orientationNumber);
            logMessage(sprintf('=> Position: %i, orientation: %i deg' ...
                ,rad2deg(positionAngle),rad2deg(orientationAngle)) ...
                ,path2LogFile,false);
            yAxis = [0 1 0];
            rotationMatrixOrientation = get3DRotationMatrix( ...
                orientationAngle,yAxis);
            totalRotationMatrix = ...
                rotationMatrixOrientation*rotationMatrixPosition;
            logMessage('    Transforming coordinates.',path2LogFile ...
                ,false);
            [rotatedX,rotatedY,rotatedZ] =  ...
                rotateTrajectoriesWithRotationMatrix( ...
                totalRotationMatrix,nearestNeighboursX ...
                ,nearestNeighboursY,nearestNeighboursZ);
            
            logMessage('    Calculation spherical coordinates.' ...
                ,path2LogFile,false);
            [polarAngle,azimuthAngle] = ... 
                transformToSphericalCoordinates(rotatedX,rotatedY ...
                ,rotatedZ);
            
            logMessage('    Calculating spherical coordinates.' ...
                ,path2LogFile,false);
            [sphericalHarmonicZerothOrder,sphericalHarmonicFirstOrder ...
                ,sphericalHarmonicSecondOrder] ...
                = calculateSphericalHarmonics(polarAngle,azimuthAngle ...
                ,nearestNeighbourDistancesPow3);
            
            logMessage('    Calculating correlation function.' ...
                ,path2LogFile,false);
            correlationFunction0W0 = calculateCorrelationFunction( ...
                sphericalHarmonicZerothOrder,lags);
            correlationFunction1W0 = calculateCorrelationFunction( ...
                sphericalHarmonicFirstOrder,lags);
            correlationFunction2W0 = calculateCorrelationFunction( ...
                sphericalHarmonicSecondOrder,lags);
            
            correlationFunction0W0Saver(orientationNumber ...
                ,positionNumber,:) = (squeeze( ...
                correlationFunction0W0Saver(orientationNumber ...
                ,positionNumber,:))'*(atomCounter-1) + ...
                correlationFunction0W0)/atomCounter;
            correlationFunction1W0Saver(orientationNumber ...
                ,positionNumber,:) = (squeeze( ...
                correlationFunction1W0Saver( ...
                orientationNumber,positionNumber,:))'*(atomCounter-1) + ...
                correlationFunction1W0)/atomCounter;
            correlationFunction2W0Saver(orientationNumber ...
                ,positionNumber,:) = (squeeze( ...
                correlationFunction2W0Saver( ...
                orientationNumber,positionNumber,:))'*(atomCounter-1) + ...
                correlationFunction2W0)/atomCounter;            
            
            logMessage('    Calculating spectral density.' ...
                ,path2LogFile,false);
            [spectralDensity1W0,spectralDensity2W0] = ...
                calculateSpectralDensities(correlationFunction1W0 ...
                ,correlationFunction2W0,omega0,deltaT,lags);
            
            logMessage('    Calculating relaxation rate.' ...
                ,path2LogFile,false);
            r1WithPerturbationTheory(orientationNumber,positionNumber ...
                ,atomCounter) = calculateR1WithSpectralDensity( ...
                spectralDensity1W0,spectralDensity2W0,dipolDipolConstant);
            logMessage(sprintf('--> R1 = %.4f.' ...
                ,r1WithPerturbationTheory(orientationNumber ...
                ,positionNumber,atomCounter)),path2LogFile);   
        end
        printDottedBreakLineToLogFile(path2LogFile);
    end
    
    if mod(atomCounter,configuration.savingInterval) == 0
        lastSavingDate = datestr(now,'yyyymmdd_HHMM');
        createDataSavingObject();
        save(path2Save,'dataSavingObject','-v7.3');
        logMessage('Saved data',path2LogFile);
    end
    logMessage(sprintf('Finished atom %i. Needed time %.4f' ...
        ,atomNumber,toc(overallForAtom)),path2LogFile,false);
    logMessage(sprintf('Calculated %i atoms.',atomCounter),path2LogFile);
    printBreakLineToLogFile(path2LogFile);
    atomCounter = atomCounter + 1;
end

