%% load configuration
clc

configuration = readConfigurationFile('config.conf');

runOnServer = configuration.runOnServer;
if runOnServer
    addpath(configuration.path2LibraryOnServer);
    path2Data = configuration.path2DataOnServer;
else
    addpath(configuration.path2LibraryOnLocalMachine)
    path2Data = configuration.path2DataOnLocalMachine;
end

%% Load Trajectories
loaded = configuration.dataLoaded;
if not(loaded)
    
    disp('Start loading data')
    
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

showFigures = configuration.showFigures;
saveFigures = configuration.saveFigures;


%% Define constants
constants = readConstantsFile(path2ConstantsFile);
picoSecond = constants.picoSecond;
hbar = constants.hbar;                              % [Js]
gammaRad = constants.gammaRad;                      % [rad/Ts]
mu0 = constants.mu0;                                % [N/A^2] Vacuum permeability
Nm = constants.Nm;                                  % [m] Nanometer (Traj. format)
DD = 3/4*(mu0/(4*pi)*hbar*gammaRad^2 )^2/(Nm^6);    % J/rad
B0 = configuration.B0;     
deltaT = configuration.deltaT;
omega0 = gammaRad*B0;                               % [rad/s]: Larmor (anglular) frequency

%% Define simulation parameters
[numberOfHs,timeSteps] = size(trajectoryX);
lags = round(configuration.fractionForLags*timeSteps);
nearestNeighbours = configuration.nearestNeighbours;
fibreOrientationsCount = configuration.fibreOrientationsCount;
positionsAtOrientationCount = configuration.myelinPositionsCount;

%% Start simulation
stopWatch = zeros(1,numberOfHs);

orientationAngles = deg2rad(linspace(0,90,fibreOrientationsCount));
positionAngles = deg2rad(linspace(0,360,positionsAtOrientationCount+1));
positionAngles = positionAngles(1:end-1);

correlationFunctionW0Saver = complex(zeros(fibreOrientationsCount ...
    ,positionsAtOrientationCount,lags));
correlationFunction2W0Saver = complex(zeros(fibreOrientationsCount ...
    ,positionsAtOrientationCount,lags));

r1WithPerturbationTheory = zeros(fibreOrientationsCount ...
    ,positionsAtOrientationCount,numberOfHs);
averageRelaxationRates = zeros(fibreOrientationsCount ...
    ,positionsAtOrientationCount,numberOfHs);
effectiveRelaxationRates = zeros(fibreOrientationsCount,numberOfHs);

meanPositions = zeros(numberOfHs,3);
atomCounter = 0;

disp('Starting Simulation.')

for atomNumber=randperm(numberOfHs)
    disp('=============================')
    disp(['Atom number: ' num2str(atomNumber)])
    disp(['Calculated ' num2str(atomCounter) ' atoms.'])
    atomCounter = atomCounter+1;
    
    meanPositions(atomNumber,:) = [mean(trajectoryX(atomNumber,:)) ...
        ,mean(trajectoryY(atomNumber,:)),mean(trajectoryZ(atomNumber,:))];
    
    [relativeXPositions,relativeYPositions,relativeZPositions ...
        ,distances] = calculatePositionsAndDistances(trajectoryX,trajectoryY ...
        ,trajectoryZ,atomNumber);
    
    [nearestNeighboursX,nearestNeighboursY,nearestNeighboursZ ...
        ,nearestNeighbourDistancesPow3] = findNearestNeighbours( ...
        distances,nearestNeighbours+1,atomNumber,relativeXPositions...
        ,relativeYPositions,relativeZPositions);
    tic;
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
        for positionNumber = 1:positionsAtOrientationCount
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
                ,atomNumber) = calculateR1WithSpectralDensity( ...
                spectralDensityW0,spectralDensity2W0,DD);
        end
        toc(positionsTic)
    end
    
    averageRelaxationRates(:,:,atomNumber) = mean( ...
        r1WithPerturbationTheory(:,:,1:atomNumber),3);
    effectiveRelaxationRates(:,atomNumber) = mean(squeeze( ...
        averageRelaxationRates(:,:,atomNumber)),2);
    
    if mod(atomNumber,10) == 0
            save(path2Save ,'r1WithPerturbationTheory' ...
                ,'averageRelaxationRates' ...
                ,'effectiveRelaxationRates' ...
                ,'correlationFunctionW0Saver' ...
                ,'correlationFunction2W0Saver' ...
                ,'meanPositions','deltaT' ...
                ,'timeSteps','lags','atomNumber','B0','-v7.3')
    end
    stopWatch(atomNumber) = toc;
    
    disp(['Atom ' num2str(atomNumber) ' done.'])
    disp(['Needed time: ' num2str(stopWatch(atomNumber)) ' seconds.'])
    
    if showFigures && ~runOnServer
        figs(1) = figure(1);
        hold on
        legendEntries = {};
        for orientationNumber = 1:size(effectiveRelaxationRates,1)
            plot(effectiveRelaxationRates(orientationNumber,1:atomNumber));
            legendEntries{orientationNumber} = num2str( ...
                rad2deg(orientationAngles(orientationNumber)), ...
                'Theta = %.2f°'); %#ok<SAGROW>
        end
        legend(legendEntries);
        title('Relaxation Rate R1')
        xlabel('Epoches')
        ylabel('R1 [Hz]')
        grid on
        drawnow
        hold off
        
        if saveFigures
           savefig(figs,path2Save); 
        end
    end
end

