%% load configuration
clc

configuration = readConfigurationFile('config.conf');

runOnServer = configuration.runOnServer;
if runOnServer
    addpath(genpath(configuration.path2LibraryOnServer));
    path2Data = configuration.path2DataOnServer;
else
    addpath(genpath(configuration.path2LibraryOnLocalMachine))
    path2Data = configuration.path2DataOnLocalMachine;
end

%% Load Trajectories
loaded = configuration.dataLoaded;
if not(loaded)
    
    disp('Loading data')
    
    fileName = configuration.fileName;
    path2File = [path2Data fileName '.mat'];
    data = load(path2File);
    hydrogenTrajectories = data.(configuration.dataFieldName);
    
    disp('Data successfully loaded')
    
    trajectoryX = squeeze(hydrogenTrajectories(:,1,:));
    trajectoryY = squeeze(hydrogenTrajectories(:,2,:));
    trajectoryZ = squeeze(hydrogenTrajectories(:,3,:));
    
end
%% Set Up
clearvars -except  trajectoryX trajectoryY trajectoryZ configuration ...
    fileName runOnServer

if runOnServer
    path2ConstantsFile = configuration.path2ConstantsFileOnServer;
    path2Save = [configuration.path2ResultsOnServer fileName ...
        configuration.resultsSuffix];
else
    path2ConstantsFile = configuration.path2ConstantsFileOnLocalMachine;
    path2Save = [configuration.path2ResultsOnLocalMachine fileName ...
        configuration.resultsSuffix];
end

%% Define constants
disp('Defining constants')
constants = readConstantsFile(path2ConstantsFile);
picoSecond = constants.picoSecond;
hbar = constants.hbar;                              % [Js]
gammaRad = constants.gammaRad;                      % [rad/Ts]
mu0 = constants.mu0;                                % [N/A^2] Vacuum permeability
Nm = constants.Nm;                                  % [m] Nanometer (Traj. format)
DD = 3/4*(mu0/(4*pi)*hbar*gammaRad^2)^2/(Nm^6);    % J/rad
B0 = configuration.B0;     
deltaT = configuration.deltaT;
omega0 = gammaRad*B0;                               % [rad/s]: Larmor (anglular) frequency

%% Define simulation parameters
disp('Defining simulation parameters')
[numberOfHs,timeSteps] = size(trajectoryX);
lags = round(configuration.fractionForLags*timeSteps);
nearestNeighbours = configuration.nearestNeighbours;
atomsToCalculate = configuration.atomsToCalculate;

%% Start simulation
stopWatch = zeros(1,atomsToCalculate);

orientationAngles = deg2rad(getValuesFromStringEnumeration( ...
    configuration.fibreOrientations,';','numeric'));
positionAngles = deg2rad(getValuesFromStringEnumeration( ...
    configuration.myelinPositions,';','numeric'));

fibreOrientationsCount = size(orientationAngles,2);
positionsInMyelinCount = size(positionAngles,2);

correlationFunctionW0Saver = complex(zeros(fibreOrientationsCount ...
    ,positionsInMyelinCount,lags));
correlationFunction2W0Saver = complex(zeros(fibreOrientationsCount ...
    ,positionsInMyelinCount,lags));

r1WithPerturbationTheory = zeros(fibreOrientationsCount ...
    ,positionsInMyelinCount,atomsToCalculate);

meanPositions = [mean(trajectoryX,2) mean(trajectoryY,2) ...
    mean(trajectoryZ,2)];

atomCounter = 0;
atomIndex = zeros(1,atomsToCalculate);

disp('Starting Simulation.')

for atomNumber=randperm(numberOfHs)
    wholeTic = tic;
    transformationTic = tic;
    disp('=============================')
    disp(['Atom number: ' num2str(atomNumber)])
    disp(['Calculated ' num2str(atomCounter) ' atoms.'])
    atomCounter = atomCounter+1;
    atomIndex(atomCounter) = atomNumber;
    
    [trajectoryX,trajectoryY,trajectoryZ,distances] ...
        = calculateRelativePositionsAndDistances(trajectoryX ...
        ,trajectoryY,trajectoryZ,atomNumber);
    
    [nearestNeighboursX,nearestNeighboursY,nearestNeighboursZ ...
        ,nearestNeighbourDistancesPow3] = findNearestNeighbours( ...
        distances,nearestNeighbours+1,atomNumber,trajectoryX ...
        ,trajectoryY,trajectoryZ);
    
    clear distances
    
    transformationTime = toc(transformationTic);
    disp(['Time for transformation of data: ' ...
        num2str(transformationTime)]);
    orientationTic = tic;
    for orientationNumber = 1:fibreOrientationsCount
        orientationAngle = orientationAngles(orientationNumber);
        yAxis = [0 1 0];
        rotationMatrixOrientation = get3DRotationMatrix( ...
            orientationAngle,yAxis);
        [rotatedX,rotatedY,rotatedZ] = ...
            rotateTrajectoriesWithRotationMatrix( ...
            rotationMatrixOrientation,nearestNeighboursX ...
            ,nearestNeighboursY,nearestNeighboursZ);
        
        positionsTic = tic;
        for positionNumber = 1:positionsInMyelinCount
            positionAngle = positionAngles(positionNumber);
            zAxis = [0 0 1];
            rotationMatrixPosition = get3DRotationMatrix(positionAngle ...
                ,zAxis);
            [rotatedX,rotatedY,rotatedZ] = ...
                rotateTrajectoriesWithRotationMatrix( ...
                rotationMatrixPosition,rotatedX,rotatedY,rotatedZ);
            
            [polarAngle,azimuthAngle] = ... 
                transformToSphericalCoordinates(rotatedX,rotatedY ...
                ,rotatedZ);
            
            [firstOrderSphericalHarmonic,secondOrderSphericalHarmonic] ...
                = calculateSphericalHarmonics(polarAngle,azimuthAngle ...
                ,nearestNeighbourDistancesPow3);
            
            correlationFunctionW0 = calculateCrossCorrelationFunction( ...
                firstOrderSphericalHarmonic ...
                ,conj(firstOrderSphericalHarmonic),lags);
            correlationFunction2W0 = calculateCrossCorrelationFunction( ...
                secondOrderSphericalHarmonic ...
                ,conj(secondOrderSphericalHarmonic),lags);
            
            correlationFunctionW0Saver(orientationNumber ...
                ,positionNumber,:) = squeeze( ...
                correlationFunctionW0Saver(orientationNumber ...
                ,positionNumber,:))'+correlationFunctionW0;
            
            correlationFunction2W0Saver(orientationNumber ...
                ,positionNumber,:) = squeeze( ...
                correlationFunction2W0Saver(orientationNumber ...
                ,positionNumber,:))'+correlationFunction2W0;
            
            [spectralDensityW0,spectralDensity2W0] = ...
                calculateSpectralDensities(correlationFunctionW0 ...
                ,correlationFunction2W0,omega0,deltaT,lags);
            
            r1WithPerturbationTheory(orientationNumber,positionNumber ...
                ,atomCounter) = calculateR1WithSpectralDensity( ...
                spectralDensityW0,spectralDensity2W0,DD);
        end
        disp(['Time for one orientation: ' num2str(toc(positionsTic))])
    end
    disp(['Time for all orientations: ' num2str(toc(orientationTic))])
    if mod(atomCounter,10) == 0
        savingTic = tic;
        save(path2Save ,'r1WithPerturbationTheory' ...
            ,'correlationFunctionW0Saver' ...
            ,'correlationFunction2W0Saver' ...
            ,'meanPositions','deltaT' ...
            ,'timeSteps','lags','atomCounter','B0' ...
            ,'orientationAngles','positionAngles' ...
            ,'atomIndex','nearestNeighbours','atomsToCalculate' ...
            ,'fileName','stopWatch'...
            ,'-v7.3')
        disp(['Time for saving variables: ' num2str(toc(savingTic))])
    end
    
    stopWatch(atomCounter) = toc(wholeTic);
    disp(['Atom ' num2str(atomNumber) ' done.'])
    disp(['Overall needed time: ' num2str(stopWatch(atomCounter)) ...
        ' seconds.'])
end

