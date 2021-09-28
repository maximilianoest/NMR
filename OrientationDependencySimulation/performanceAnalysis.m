%% Set Up System
clc
clear all %#ok<CLALL>

setUpTimer = tic;
configuration = readConfigurationFile('config.txt');
if configuration.runOnServer
    addpath(genpath(configuration.path2LibraryOnServer));
else
    addpath(genpath(configuration.path2LibraryOnLocalMachine))
end
timeTracks = {};

[path2Data,path2Save,path2ConstantsFile,path2LogFile] = ...
    setUpSystemBasedOnMachine(configuration);
deleteLogFile(path2LogFile);
fileName = configuration.fileName;

logMessage('System is set up.', path2LogFile);
logMemoryUsage(path2LogFile);

%% Load data
tic;
logMessage('Start loading data.', path2LogFile);
[trajectoryX,trajectoryY,trajectoryZ] = loadTrajectoriesFromData( ...
    configuration,path2Data);
timeTracks.dataLoad = toc;
logMessage(sprintf('Data successfully loaded in %.3f s.', ...
    timeTracks.dataLoad), path2LogFile);
logMemoryUsage(path2LogFile);

%% Define constants
logMessage('Defining constants.',path2LogFile,false);
deltaT = configuration.deltaT;
constants = readConstantsFile(path2ConstantsFile);

dipolDipolConstant = 3/4*(constants.vaccumPermeability/(4*pi) ...
    *constants.hbar*constants.gyromagneticRatioOfHydrogenAtom^2)^2 ...
    /(constants.nanoMeter^6);
omega0 = constants.gyromagneticRatioOfHydrogenAtom ...
    *configuration.mainMagneticField; 
logMemoryUsage(path2LogFile);

%% Define simulation parameters
logMessage('Defining simulation parameters.',path2LogFile,false);

[numberOfHs,timeSteps] = size(trajectoryX);
logMessage(sprintf(['    Found %d hydrogen atoms at %d time steps of ' ...
    '%.3d s'],numberOfHs,timeSteps,deltaT),path2LogFile,false);

lags = round(configuration.fractionForLags*timeSteps);
logMessage(sprintf(['    The lag is set to %d time steps, resulting ' ...
    'in a simulation time of %d s. NO FUNCTIONALITY'], lags, lags*deltaT),path2LogFile ...
    ,false);

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

shiftForCorrelationFunction = configuration.shiftForCorrelationFunction;
tmp = false(1,lags);
tmp = tmp(1:shiftForCorrelationFunction:end);
correlationFunction1W0Saver = zeros(fibreOrientationsCount ...
    ,positionsInMyelinCount,atomsToCalculate,length(tmp) ...
    ,'like',single(1j));
correlationFunction2W0Saver = zeros(fibreOrientationsCount ...
    ,positionsInMyelinCount,atomsToCalculate,length(tmp) ...
    ,'like',single(1j));
clearvars tmp
logMessage('    Created correlation function saver.' ...
    , path2LogFile,false);
logMemoryUsage(path2LogFile);

%% Start simulation
logMessage('Preallocation of some other arrays.',path2LogFile,false);

r1WithPerturbationTheory = zeros(fibreOrientationsCount ...
    ,positionsInMyelinCount,atomsToCalculate);
meanPositions = single([mean(trajectoryX,2) mean(trajectoryY,2) ...
    mean(trajectoryZ,2)]);

atomCounter = 1;
atomIndex = zeros(1,atomsToCalculate);
timeTracks.setUp = toc(setUpTimer);
calculationSteps = fibreOrientationsCount*positionsInMyelinCount;

randomSequenceOfAtoms = randperm(numberOfHs);

relativeX = zeros(numberOfHs,timeSteps);
relativeY= zeros(numberOfHs,timeSteps);
relativeZ = zeros(numberOfHs,timeSteps);

nearestNeighboursX = zeros(nearestNeighbours,timeSteps);
nearestNeighboursY = zeros(nearestNeighbours,timeSteps);
nearestNeighboursZ = zeros(nearestNeighbours,timeSteps);
nearestNeighbourDistancesPow3 = zeros(nearestNeighbours,timeSteps);

rotatedX = zeros(nearestNeighbours,timeSteps);
rotatedY = zeros(nearestNeighbours,timeSteps);
rotatedZ = zeros(nearestNeighbours,timeSteps);

polarAngle = zeros(nearestNeighbours,timeSteps);
azimuthAngle = zeros(nearestNeighbours,timeSteps);

firstOrderSphericalHarmonic = zeros(nearestNeighbours,timeSteps);
secondOrderSphericalHarmonic = zeros(nearestNeighbours,timeSteps);

correlationFunction1W0 = complex(zeros(1,lags));
correlationFunction2W0 = complex(zeros(1,lags));

timeTracks.relativePositions = zeros(1,atomsToCalculate);
timeTracks.nearestNeighbours = zeros(1,atomsToCalculate);
timeTracks.transformation = zeros(1,atomsToCalculate);
timeTracks.overallForAtom = zeros(1,atomsToCalculate);
timeTracks.transformation = zeros(1,atomsToCalculate);
timeTracks.sphericalCoordinates = zeros(1,atomsToCalculate);
timeTracks.sphericalHarmonics = zeros(1,atomsToCalculate);
timeTracks.correlationFunction = zeros(1,atomsToCalculate);
timeTracks.spectralDensity = zeros(1,atomsToCalculate);
timeTracks.relaxationRate = zeros(1,atomsToCalculate);
timeTracks.position = zeros(1,atomsToCalculate);

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

logMessage(sprintf('Starting simulation after %f s for set up.' ...
    ,timeTracks.setUp),path2LogFile);
logMemoryUsage(path2LogFile);
printBreakLineToLogFile(path2LogFile);



% for atomNumber = randomSequenceOfAtoms(atomCounter:end)
for atomNumber = 1:numberOfHs
    overallForAtom = tic;
    logMessage(sprintf('Selected atom number %i',atomNumber),path2LogFile);
    atomIndex(atomCounter) = atomNumber;
    
    relativePositionsTimer = tic;
    [relativeX,relativeY,relativeZ] ...
        = calculateRelativePositions(trajectoryX ...
        ,trajectoryY,trajectoryZ,atomNumber);
    timeTracks.relativePositions(atomCounter) = ...
        toc(relativePositionsTimer);
    logMessage('    Relative positions.',path2LogFile,false);
    logMemoryUsage(path2LogFile);
    
    nearestNeighboursTimer = tic;
    [nearestNeighboursX,nearestNeighboursY,nearestNeighboursZ ...
        ,nearestNeighbourDistancesPow3] = findNearestNeighbours( ...
        nearestNeighbours+1,atomNumber,relativeX,relativeY ...
        ,relativeZ);
    timeTracks.nearestNeighbours(atomCounter) = ...
        toc(nearestNeighboursTimer);
    logMessage('    Nearest neighbours.',path2LogFile,false);
    logMemoryUsage(path2LogFile);
    
    for positionNumber = 1:positionsInMyelinCount
        positionTimer = tic;
        positionAngle = positionAngles(positionNumber);
        zAxis = [0 0 1];
        rotationMatrixPosition = get3DRotationMatrix( ...
            positionAngle,zAxis);
        for orientationNumber = 1:fibreOrientationsCount
            transformationTimer = tic;
            orientationAngle = orientationAngles(orientationNumber);
            logMessage(sprintf('=> Orientation: %i°, position: %i°' ...
                ,rad2deg(orientationAngle),rad2deg(positionAngle)) ...
                ,path2LogFile);
            yAxis = [0 1 0];
            rotationMatrixOrientation = get3DRotationMatrix( ...
                orientationAngle,yAxis);
            totalRotationMatrix = ...
                rotationMatrixOrientation*rotationMatrixPosition;
            [rotatedX,rotatedY,rotatedZ] =  ...
                rotateTrajectoriesWithRotationMatrix( ...
                totalRotationMatrix,nearestNeighboursX ...
                ,nearestNeighboursY,nearestNeighboursZ);
            timeTracks.transformation(atomCounter) = ...
                timeTracks.transformation(atomCounter) ...
                + toc(transformationTimer);
            logMessage('    Coordinate transformation.',path2LogFile ...
                ,false);
            logMemoryUsage(path2LogFile);
            
            sphericalCoordinatesTimer = tic;
            [polarAngle,azimuthAngle] = ... 
                transformToSphericalCoordinates(rotatedX,rotatedY ...
                ,rotatedZ);
            timeTracks.sphericalCoordinates(atomCounter) = ...
                timeTracks.sphericalCoordinates(atomCounter) ...
                + toc(sphericalCoordinatesTimer);
            logMessage('    Spherical coordinates.',path2LogFile,false);
            logMemoryUsage(path2LogFile);
            
            sphericalHarmonicsTimer = tic;
            [firstOrderSphericalHarmonic,secondOrderSphericalHarmonic] ...
                = calculateSphericalHarmonics(polarAngle,azimuthAngle ...
                ,nearestNeighbourDistancesPow3);
            timeTracks.sphericalHarmonics(atomCounter) = ...
                timeTracks.sphericalHarmonics(atomCounter) ...
                + toc(sphericalHarmonicsTimer);
            logMessage('    Spherical harmonics.',path2LogFile,false);
            logMemoryUsage(path2LogFile);
            
            correlationFunctionTimer = tic;
            correlationFunction1W0 = calculateCorrelationFunction( ...
                firstOrderSphericalHarmonic,lags);
            correlationFunction2W0 = calculateCorrelationFunction( ...
                secondOrderSphericalHarmonic,lags);
            correlationFunction1W0Saver(orientationNumber ...
                ,positionNumber,atomCounter,:) = ...
                correlationFunction1W0(1:shiftForCorrelationFunction:end);
            correlationFunction2W0Saver(orientationNumber ...
                ,positionNumber,atomCounter,:) = ...
                correlationFunction2W0(1:shiftForCorrelationFunction:end);
            timeTracks.correlationFunction(atomCounter) = ...
                timeTracks.correlationFunction(atomCounter) ...
                + toc(correlationFunctionTimer);
            logMessage('    Correlation function.',path2LogFile,false);
            logMemoryUsage(path2LogFile);
            
            spectralDensityTimer = tic;
            [spectralDensityW0,spectralDensity2W0] = ...
                calculateSpectralDensities(correlationFunction1W0 ...
                ,correlationFunction2W0,omega0,deltaT,lags);
            timeTracks.spectralDensity(atomCounter) = ...
                timeTracks.spectralDensity(atomCounter) ...
                + toc(spectralDensityTimer);
            logMessage('    Spectral density.',path2LogFile,false);
            logMemoryUsage(path2LogFile);
            
            relaxationRateTimer = tic;
            r1WithPerturbationTheory(orientationNumber,positionNumber ...
                ,atomCounter) = calculateR1WithSpectralDensity( ...
                spectralDensityW0,spectralDensity2W0,dipolDipolConstant);
            timeTracks.relaxationRate(atomCounter) = ...
                timeTracks.relaxationRate(atomCounter) ...
                + toc(relaxationRateTimer);
            logMessage('    Relaxation rate.',path2LogFile,false);
            logMemoryUsage(path2LogFile);
            logMessage(sprintf('---> R1 = %.4f.' ...
                ,r1WithPerturbationTheory(orientationNumber ...
                ,positionNumber,atomCounter)),path2LogFile);   
        end
        timeTracks.position(atomCounter) = ...
            timeTracks.position(atomCounter) + toc(positionTimer);
        logMessage(sprintf(['Finished position %i°. Needed '...
            'time %.4f'],rad2deg(positionAngle),toc(positionTimer)) ...
            ,path2LogFile,false);
    end
    logMemoryUsage(path2LogFile);
    timeTracks.overallForAtom(atomCounter) = toc(overallForAtom);
    timeTracks.transformation(atomCounter) = ...
        timeTracks.transformation(atomCounter)/calculationSteps;
    timeTracks.sphericalCoordinates(atomCounter) = ...
        timeTracks.sphericalCoordinates(atomCounter)/calculationSteps;
    timeTracks.sphericalHarmonics(atomCounter) = ...
        timeTracks.sphericalHarmonics(atomCounter)/calculationSteps;
    timeTracks.correlationFunction(atomCounter) = ...
        timeTracks.correlationFunction(atomCounter)/calculationSteps;
    timeTracks.spectralDensity(atomCounter) = ...
        timeTracks.spectralDensity(atomCounter)/calculationSteps;
    timeTracks.relaxationRate(atomCounter) = ...
        timeTracks.relaxationRate(atomCounter)/calculationSteps;
    timeTracks.position(atomCounter) = ...
        timeTracks.position(atomCounter)/positionsInMyelinCount;
    
    if mod(atomCounter,1) == 0
        lastSavingDate = datestr(now,'yyyymmdd_HHMM');
        save(path2Save ,'path2Data','path2Save','path2ConstantsFile' ...
            ,'path2LogFile','timeTracks'...
            ,'r1WithPerturbationTheory','configuration' ...
            ,'correlationFunction1W0Saver' ...
            ,'correlationFunction2W0Saver','meanPositions','deltaT' ...
            ,'timeSteps','lags','atomCounter','orientationAngles' ...
            ,'positionAngles','atomsToCalculate' ...
            ,'atomIndex','nearestNeighbours','atomsToCalculate' ...
            ,'fileName','numberOfHs' ...
            ,'startDateOfSimulation','lastSavingDate','constants' ...
            ,'shiftForCorrelationFunction','dipolDipolConstant' ...
            ,'constants','lags','randomSequenceOfAtoms','-v7.3')
        logMessage('Saved data',path2LogFile);
    end
    logPerfomanceData(atomNumber,atomCounter,path2LogFile,timeTracks)
    logMessage(sprintf('Calculated %i atoms',atomCounter),path2LogFile);
    printDottedBreakLineToLogFile(path2LogFile);
    atomCounter = atomCounter + 1;
end

