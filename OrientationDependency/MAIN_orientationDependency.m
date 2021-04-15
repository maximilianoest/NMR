%% load configuration

configuration = readConfigurationFile('config.conf');

runOnServer = configuration.runOnServer;
if runOnServer
    addpath(configuration.path2LibraryOnServer);
else
    addpath(configuration.path2LibraryOnLocalMachine)
end

%% Load Trajectories
loaded = configuration.dataLoaded;
if not(loaded)
    
    clc
    clearvars -except configuration
    close all
    
    LIPID = configuration.lipid;
    path2Data = configuration.path2Data;
    
    switch LIPID
        case{'DOPS'}
            fileName = configuration.fileNameDOPS;
        case{'PLPC'}
            fileName = configuration.fileNamePLPC;
        case{'PSM'}
            fileName = configuration.fileNamePSM;
        case{'Myelin'}
            fileName = configuration.fileNameMyelin;
        otherwise
            warning(['The lipid "' LIPID '" does not exist'])
    end
    
    path2File = [path2Data fileName '.mat'];
    data = load(path2File);
    hydrogenTrajectories = data.(configuration.dataFieldName);
    
    disp('Data successfully loaded')
    
    trajX = squeeze(hydrogenTrajectories(:,1,:));
    trajY = squeeze(hydrogenTrajectories(:,2,:));
    trajZ = squeeze(hydrogenTrajectories(:,3,:));
    
end
%% Set Up
clearvars -except  trajX trajY trajZ configuration fileName LIPID

calculateSchroedingerEquation = configuration.calculateSchroedingerEquation;
path2ConstantsFile = configuration.path2ConstantsFile;
path2Save = [configuration.path2Results fileName ...
    configuration.resultsSuffix];
showFigures = configuration.showFigures;
outputLogFileName = configuration.outputLogFileName;
runOnServer = configuration.runOnServer;


%% Define constants
constants = readConstantsFile(path2ConstantsFile);
picoSecond = constants.picoSecond;
deltaT = configuration.deltaT*picoSecond;
hbar = constants.hbar;                              % [Js]
gammaRad = constants.gammaRad;                      % [rad/Ts]
B0 = constants.B0;                                  % [T]
mu0 = constants.mu0;                                % [N/A^2] Vacuum permeability
Nm = constants.Nm;                                  % [m] Nanometer (Traj. format)
DD = 3/4*(mu0/(4*pi)*hbar*gammaRad^2 )^2/(Nm^6);    % J/rad
omega0 = gammaRad*B0;                               % [rad/s]: Larmor (anglular) frequency

%% Define simulation parameters
[numberOfHs,timeSteps] = size(trajX);
lags = round(configuration.fractionForLags*timeSteps);
nearestNeighbours = configuration.nearestNeighbours;
fibreOrientationsCount = configuration.fibreOrientationsCount;
positionsAtOrientationCount = configuration.myelinPositionsCount;

%% Start simulation
orientationAngles = deg2rad(linspace(0,90,fibreOrientationsCount));
positionAngles = deg2rad(linspace(0,360,positionsAtOrientationCount+1));
positionAngles = positionAngles(1:end-1);

correlationFunctionW0Saver = [];
correlationFunction2W0Saver = [];
r1WithPerturbationTheory = zeros(fibreOrientationsCount ...
    ,positionsAtOrientationCount,numberOfHs);

averageRelaxationRates = zeros(fibreOrientationsCount ...
    ,positionsAtOrientationCount,numberOfHs);
effectiveRelaxationRates = zeros(fibreOrientationsCount,numberOfHs);

meanPositions = zeros(numberOfHs,3);
molecules=1:numberOfHs;
%molecules=molecules(randperm(length(molecules)));

if showFigures
    figure('Name','Relaxation Rates')
end

disp('Start Simulation.')

for atomNumber=1:numberOfHs
    if runOnServer
        disp(['Atom number: ' num2str(atomNumber)])
        disp('=============================')
    else
        fileId = fopen(outputLogFileName,'a');
        fprintf(fileId,['========================== \n' ...
            ,'  Atom Number: %d \n'],atomNumber);
        fclose(fileId);
    end
    
    meanPositions(atomNumber,:) = [mean(trajX(atomNumber,:)) ...
        ,mean(trajY(atomNumber,:)),mean(trajZ(atomNumber,:))];
    
    [relativeXPositions,relativeYPositions,relativeZPositions ...
        ,distances] = calculatePositionsAndDistances(trajX,trajY ...
        ,trajZ,atomNumber);
    
    [nearestNeighboursX,nearestNeighboursY,nearestNeighboursZ ...
        ,nearestNeighbourDistancesPow3] = findNearestNeighbours( ...
        distances,nearestNeighbours+1,atomNumber,relativeXPositions...
        ,relativeYPositions,relativeZPositions);
    
    for orientationNumber = 1:fibreOrientationsCount
        orientationAngle = orientationAngles(orientationNumber);
        if runOnServer
            disp('----------------------------')
            disp(['Orientation: ' num2str(rad2deg(orientationAngle))])
        else
            fileId = fopen(outputLogFileName,'a');
            fprintf(fileId,['--------------------- \n' ...
                'Orientation: %f \n'],rad2deg(orientationAngle));
            fclose(fileId);
        end
        yAxis = [0 1 0];
        rotationMatrixOrientation = get3DRotationMatrix( ...
            orientationAngle,yAxis);
        [rotatedX,rotatedY,rotatedZ] = ...
            rotateTrajectoriesWithRotationMatrix( ...
            rotationMatrixOrientation,nearestNeighboursX ...
            ,nearestNeighboursY,nearestNeighboursZ);
        
        for positionNumber = 1:positionsAtOrientationCount
            positionAngle = positionAngles(positionNumber);
            if runOnServer
                disp(['Position: ' num2str(rad2deg(positionAngle))])
            else
                fileId = fopen(outputLogFileName,'a');
                fprintf(fileId,'Position: %f \n',rad2deg(positionAngle));
                fclose(fileId);
            end
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
                ,positionNumber,atomNumber,:) = sum(correlationFunctionW0);  
            correlationFunction2W0Saver(orientationNumber ...
                ,positionNumber,atomNumber,:) = ...
                sum(correlationFunction2W0);
            
            [spectralDensityW0,spectralDensity2W0] = ...
                calculateSpectralDensities(correlationFunctionW0 ...
                ,correlationFunction2W0,omega0,deltaT,lags);
            
            r1WithPerturbationTheory(orientationNumber,positionNumber ...
                ,atomNumber) = calculateR1WithSpectralDensity( ...
                spectralDensityW0,spectralDensity2W0,DD);
        end
    end
    
    averageRelaxationRates(:,:,atomNumber) = mean( ...
        r1WithPerturbationTheory(:,:,1:atomNumber),3);
    effectiveRelaxationRates(:,atomNumber) = mean(squeeze( ...
        averageRelaxationRates(:,:,atomNumber)),2);
    
    if mod(atomNumber,2) == 0
            save(path2Save ,'r1WithPerturbationTheory' ...
                ,'averageRelaxationRates' ...
                ,'effectiveRelaxationRates' ...
                ,'correlationFunctionW0Saver' ...
                ,'correlationFunction2W0Saver' ...
                ,'meanPositions','deltaT' ...
                ,'timeSteps','lags','atomNumber','B0','-v7.3')
    end
    
    % FIXME: change plots to fitting relaxation rates
    if showFigures && ~runOnServer
        plot(r1WithPerturbationTheory(1:atomNumber),'b','LineWidth',1.5)
        hold on
        plot(meanR1WithPerturbationTheory(1:atomNumber),'--b' ...
            ,'LineWidth',1.5)
        plot(r1WithLipariSzabo(1:atomNumber),'k','LineWidth',1.5)
        plot(meanR1WithLipariSzabo(1:atomNumber),'--k' ...
            ,'LineWidth',1.5)
        if calculateSchroedingerEquation
            plot(r1WithSchroedingerEquation(1:atomNumber),'r' ...
                ,'LineWidth',1.5)
            plot(meanR1WithSchroedingerEquation(1:atomNumber) ...
                ,'--r','LineWidth',1.5)
            legend('Spectral density','Mean Spectral Density' ...
                ,'Lipari Szabo','Mean Lipari Szabo' ...
                ,'Schroedinger Equation','Mean Schroedinger Equ.');
        else
            legend('Spectral density','Mean Spectral Density' ...
                ,'Lipari Szabo','Mean Lipari Szabo');
        end
        title('Relaxation Rate R1')
        xlabel('Epoches')
        ylabel('R1 [Hz]')
        grid on
        drawnow
        hold off
    end
end

