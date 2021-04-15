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
    
    disp('Data loaded')
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

% orientations = [];
% positions = [];
% for angleNumber = 1:length(positionAngles)
%     orientations = [orientations orientationAngles]; %#ok<AGROW>
%     positions = [positions ones(1,length(orientationAngles)) ...
%         *positionAngles(angleNumber)]; %#ok<AGROW>
% end

correlationFunctionW0Saver = [];
correlationFunction2W0Saver = [];
r1WithPerturbationTheory = zeros(fibreOrientationsCount ...
    ,positionsAtOrientationCount,numberOfHs);


% r1WithSchroedingerEquation = zeros(1,numberOfHs);
% r1WithLipariSzabo = zeros(1,numberOfHs);

averageRelaxationRates = zeros(fibreOrientationsCount ...
    ,positionsAtOrientationCount,numberOfHs);
effectiveRelaxationRates = zeros(fibreOrientationsCount,numberOfHs);

% meanR1WithSchroedingerEquation = zeros(1,numberOfHs);
% meanR1WithLipariSzabo = zeros(1,numberOfHs);

meanPositions = zeros(numberOfHs,3);
molecules=1:numberOfHs;
%molecules=molecules(randperm(length(molecules)));

if showFigures
    figure('Name','Relaxation Rates')
end

disp('start simulation')
for atomNumber=1:numberOfHs
    fileId = fopen(outputLogFileName,'a');
    fprintf(fileId,['========================== \n' ...
        ,'  Atom Number: %d \n'],atomNumber);
    fclose(fileId);
    meanPositions(atomNumber,:) = [mean(trajX(atomNumber,:)) ...
        ,mean(trajY(atomNumber,:)),mean(trajZ(atomNumber,:))];
    
    [relativeXPositions,relativeYPositions,relativeZPositions ...
        ,distances] = calculatePositionsAndDistances(trajX,trajY ...
        ,trajZ,atomNumber);
    
    [nearestNeighboursX,nearestNeighboursY,nearestNeighboursZ ...
        ,nearestNeighbourDistancesPow3] = findNearestNeighbours( ...
        distances,nearestNeighbours+1,atomNumber,relativeXPositions...
        ,relativeYPositions,relativeZPositions);
    
    parfor orientationNumber = 1:fibreOrientationsCount
        orientationAngle = orientationAngles(orientationNumber);
        yAxis = [0 1 0];
        rotationMatrixOrientation = get3DRotationMatrix( ...
            orientationAngle,yAxis);
        [rotatedX,rotatedY,rotatedZ] = ...
            rotateTrajectoriesWithRotationMatrix( ...
            rotationMatrixOrientation,nearestNeighboursX ...
            ,nearestNeighboursY,nearestNeighboursZ);
        
        for positionNumber = 1:positionsAtOrientationCount
            positionAngle = positionAngles(positionNumber); %#ok<PFBNS>
            zAxis = [0 0 1];
            rotationMatrixPosition = get3DRotationMatrix(positionAngle ...
                ,zAxis);
            [rotatedX,rotatedY,rotatedZ] = ...
                rotateTrajectoriesWithRotationMatrix( ...
                rotationMatrixPosition,rotatedX,rotatedY,rotatedZ);
            
            
            [polarAngle,azimuthAngle] = transformToSphericalCoordinates( ...
                nearestNeighboursX,nearestNeighboursY,nearestNeighboursZ);
            
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
    
    if showFigures
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

 %     parfor orientationAndPosition = 1:length(orientations)
    %         positionAngle = positions(orientationAndPosition);
    %         zAxis = [0 0 1];
    %         rotationMatrixForPosition = get3DRotationMatrix(positionAngle ...
    %             ,zAxis);
    %         [rotatedX,rotatedY,rotatedZ] = ...
    %             rotateTrajectoriesWithRotationMatrix(rotationMatrixPosition ...
    %             ,nearestNeighboursX,nearestNeighboursY,nearestNeighboursZ);
    %
    %         orientationAngle = orientations(orientationAndPosition);
    %         yAxis = [0 1 0];
    %         rotationMatrixOrientation = get3DRotationMatrix(orientationAngle ...
    %             ,yAxis);
    %         [rotatedX,rotatedY,rotatedZ] = ...
    %             rotateTrajectoriesWithRotationMatrix( ...
    %             rotationMatrixOrientation,rotatedX,rotatedY,rotatedZ);
    %
    %         [polarAngle,azimutalAngle] = transformToSphericalCoordinates( ...
    %             nearestNeighboursX,nearestNeighboursY,nearestNeighboursZ);
    %
    %         [firstOrderSphericalHarmonic,secondOrderSphericalHarmonic] ...
    %             = calculateSphericalHarmonics(polarAngle,azimutalAngle ...
    %             ,nearestNeighbourDistancesPow3);
    %
    %         correlationFunctionW0 = calculateCrossCorrelationFunction( ...
    %             firstOrderSphericalHarmonic ...
    %             ,conj(firstOrderSphericalHarmonic),lags-1);
    %         correlationFunction2W0 = calculateCrossCorrelationFunction( ...
    %             secondOrderSphericalHarmonic ...
    %             ,conj(secondOrderSphericalHarmonic),lags-1);
    %
    %         [spectralDensityW0,spectralDensity2W0] = ...
    %             calculateSpectralDensities(correlationFunctionW0 ...
    %             ,correlationFunction2W0,omega0,deltaT,lags);
    %
    %         r1WithPerturbationTheory(atomNumber) = ...
    %             calculateR1WithSpectralDensity(spectralDensityW0 ...
    %             ,spectralDensity2W0,DD);
    %         meanR1WithPerturbationTheory(atomNumber) = mean( ...
    %             r1WithPerturbationTheory(1:atomNumber));
    %
    %     end
    
    
    
    %     %     spectralDensityLiSz1 = estimateSpectralDensityWithLipariSzaboFit( ...
%     %         correlationFunction1,1*omega0,deltaT,outputLogFileName,LIPID);
%     %     spectralDensityLiSz2 = estimateSpectralDensityWithLipariSzaboFit( ...
%     %         correlationFunction2,2*omega0,deltaT,outputLogFileName,LIPID);
%     
%     
%     
%     r1WithLipariSzabo(atomNumber) = ...
%         calculateR1WithSpectralDensity(spectralDensityLiSz1 ...
%         ,spectralDensityLiSz2,DD);
%     meanR1WithLipariSzabo(atomNumber) = ...
%         mean(r1WithLipariSzabo(1:atomNumber));
%     
%     
%     
%     if calculateSchroedingerEquation
%         r1WithSchroedingerEquation(atomNumber) = ...
%             calculateR1WithSchroedingerEquation(polarAngle,azimuthAngle ...
%             ,nearestNeighbourDistancesPow3,deltaT,path2ConstantsFile);
%         meanR1WithSchroedingerEquation(atomNumber) = ...
%             mean(r1WithSchroedingerEquation(1:atomNumber));
%     end
