%% Set Up System
clc

configuration = readConfigurationFile('config.conf');

if configuration.runOnServer
    addpath(genpath(configuration.path2LibraryOnServer));
else
    addpath(genpath(configuration.path2LibraryOnLocalMachine))
end

[path2Data,path2Save,path2ConstantsFile] = setUpSystemBasedOnMachine( ...
    configuration);
[trajectoryX,trajectoryY,trajectoryZ] = loadTrajectoriesFromData( ...
    configuration,path2Data);
fileName = configuration.fileName;

clearvars -except  trajectoryX trajectoryY trajectoryZ configuration ...
    fileName path2Save path2ConstantsFile

%% Define constants
disp('Defining constants')
deltaT = configuration.deltaT;
constants = readConstantsFile(path2ConstantsFile);

picoSecond = constants.picoSecond;
hbar = constants.hbar;
gammaRad = constants.gammaRad;
B0 = configuration.B0;
mu0 = constants.mu0;
Nm = constants.Nm;
DD = 3/4*(mu0/(4*pi)*hbar*gammaRad^2)^2/(Nm^6);
omega0 = gammaRad*B0; 

%% Define simulation parameters
disp('Defining simulation parameters')
[numberOfHs,timeSteps] = size(trajectoryX);
lags = round(configuration.fractionForLags*timeSteps);
nearestNeighbours = configuration.nearestNeighbours;
atomsToCalculate = configuration.atomsToCalculate;
startDateOfSimulation = datestr(now,'yyyymmdd');

%% Start simulation
stopWatch = zeros(1,atomsToCalculate);

orientationAngles = deg2rad(getValuesFromStringEnumeration( ...
    configuration.fibreOrientations,';','numeric'));
positionAngles = deg2rad(getValuesFromStringEnumeration( ...
    configuration.myelinPositions,';','numeric'));

fibreOrientationsCount = size(orientationAngles,2);
positionsInMyelinCount = size(positionAngles,2);

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

r1WithPerturbationTheory = zeros(fibreOrientationsCount ...
    ,positionsInMyelinCount,atomsToCalculate);

meanPositions = single([mean(trajectoryX,2) mean(trajectoryY,2) ...
    mean(trajectoryZ,2)]);

atomCounter = 0;
atomIndex = zeros(1,atomsToCalculate);

disp('Starting Simulation.')

for atomNumber = randperm(numberOfHs)
    wholeTic = tic;
    transformationTic = tic;
    disp('=============================')
    disp(['Atom number: ' num2str(atomNumber)])
    disp(['Calculated ' num2str(atomCounter) ' atoms.'])
    atomCounter = atomCounter+1;
    atomIndex(atomCounter) = atomNumber;
    
    [relativeX,relativeY,relativeZ] ...
        = calculateRelativePositions(trajectoryX ...
        ,trajectoryY,trajectoryZ,atomNumber);
    
    [nearestNeighboursX,nearestNeighboursY,nearestNeighboursZ ...
        ,nearestNeighbourDistancesPow3] = findNearestNeighbours( ...
        nearestNeighbours+1,atomNumber,relativeX,relativeY ...
        ,relativeZ);
    
    clear relativeX relativeY relativeZ
    
    transformationTime = toc(transformationTic);
    disp(['Time needed for transforming to relative variables: ' ...
        num2str(transformationTime) 's']);
    positionsTic = tic;
    for positionNumber = 1:positionsInMyelinCount
        positionAngle = positionAngles(positionNumber);
        zAxis = [0 0 1];
        rotationMatrixPosition = get3DRotationMatrix( ...
            positionAngle,zAxis);
        orientationsTic = tic;
        for orientationNumber = 1:fibreOrientationsCount
            orientationAngle = orientationAngles(orientationNumber);
            yAxis = [0 1 0];
            rotationMatrixOrientation = get3DRotationMatrix( ...
                orientationAngle,yAxis);
            totalRotationMatrix = ...
                rotationMatrixOrientation*rotationMatrixPosition;
            [rotatedX,rotatedY,rotatedZ] =  ...
                rotateTrajectoriesWithRotationMatrix( ...
                totalRotationMatrix,nearestNeighboursX ...
                ,nearestNeighboursY,nearestNeighboursZ);
            
            [r1WithPerturbationTheory( ...
                    orientationNumber,positionNumber,atomCounter) ...
                ,correlationFunction1W0Saver(orientationNumber ...
                    ,positionNumber,atomCounter,:) ...
                ,correlationFunction2W0Saver(orientationNumber ...
                    ,positionNumber,atomCounter,:)] = ...
                calculateR1RatesAndCorrelationFunctions(rotatedX ...
                ,rotatedY,rotatedZ,nearestNeighbourDistancesPow3 ...
                ,lags,omega0,deltaT,DD,shiftForCorrelationFunction);
            
        end
        disp(['Time needed for all orientations at position ' ...
            num2str(rad2deg(positionAngle)) ': ' ...
            num2str(toc(orientationsTic)) 's'])
    end
    disp(['Time needed for all variations : ' ...
        num2str(toc(positionsTic)) 's'])
    if mod(atomCounter,5) == 0
        disp('Start saving data.')
        savingTic = tic;
        lastSavingDate = datestr(now,'yyyymmdd_HHMM');
        save(path2Save ,'r1WithPerturbationTheory' ...
            ,'correlationFunction1W0Saver' ...
            ,'correlationFunction2W0Saver' ...
            ,'meanPositions','deltaT' ...
            ,'timeSteps','lags','atomCounter','B0' ...
            ,'orientationAngles','positionAngles' ...
            ,'atomIndex','nearestNeighbours','atomsToCalculate' ...
            ,'fileName','stopWatch','numberOfHs' ...
            ,'startDateOfSimulation','lastSavingDate','constants' ...
            ,'shiftForCorrelationFunction','DD','-v7.3')
        disp(['Time needed for saving variables: ' ...
            num2str(toc(savingTic)) 's'])
    end
    
    stopWatch(atomCounter) = toc(wholeTic);
    disp(['Atom ' num2str(atomNumber) ' done.'])
    disp(['Overall time needed for one Atom: ' ...
        num2str(stopWatch(atomCounter)) 's'])
end

