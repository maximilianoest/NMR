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
        otherwise
            warning(['The lipid "' LIPID '" does not exist'])
    end
    
    path2File = [path2Data fileName];
    loadReducedDataSet = configuration.loadReducedDataSet;
    if loadReducedDataSet
         path2File = [path2File '_reduced.mat'];
    else
         path2File = [path2File '.mat' ];
    end
    
    data = load(path2File);

    if loadReducedDataSet
        data = data(1).DataRed;
    end
     
    lipidHydrogens = data.(fileName);

    trajX = squeeze(lipidHydrogens(:,1,:));
    trajY = squeeze(lipidHydrogens(:,2,:));
    trajZ = squeeze(lipidHydrogens(:,3,:));
end
%% Set Up
clearvars -except  trajX trajY trajZ configuration fileName LIPID

calculateSchroedingerEquation = configuration.calculateSchroedingerEquation;
path2ConstantsFile = configuration.path2ConstantsFile;
path2Save = [configuration.path2Results fileName '_resultsRelaxationRates'];
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

 %% Start simulation

r1WithPerturbationTheory = zeros(1,numberOfHs);
r1WithSchroedingerEquation = zeros(1,numberOfHs);
r1WithLipariSzabo = zeros(1,numberOfHs);

meanR1WithPerturbationTheory = zeros(1,numberOfHs);
meanR1WithSchroedingerEquation = zeros(1,numberOfHs);
meanR1WithLipariSzabo = zeros(1,numberOfHs);

meanPositions = zeros(numberOfHs,3);
molecules=1:numberOfHs;
%molecules=molecules(randperm(length(molecules)));

if showFigures
    figure('Name','Relaxation Rates')
end

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
    
    [theta,phi] = transformToSphericalCoordinates(nearestNeighboursX ...
        ,nearestNeighboursY,nearestNeighboursZ);
    
    [F1,F2] = calculateSphericalHarmonics(theta,phi ...
        ,nearestNeighbourDistancesPow3);
    
    correlationFunction1 = calculateCrossCorrelationFunction(F1 ...
        ,conj(F1),lags-1);
    correlationFunction2 = calculateCrossCorrelationFunction(F2 ...
        ,conj(F2),lags-1);
    
    spectralDensityLiSz1 = estimateSpectralDensityWithLipariSzaboFit( ...
        correlationFunction1,1*omega0,deltaT,outputLogFileName,LIPID);
    spectralDensityLiSz2 = estimateSpectralDensityWithLipariSzaboFit( ...
        correlationFunction2,2*omega0,deltaT,outputLogFileName,LIPID);
    
    [spectralDensity1,spectralDensity2] = calculateSpectralDensities( ...
        correlationFunction1,correlationFunction2,omega0,deltaT,lags);
    
    r1WithLipariSzabo(atomNumber) = ...
        calculateR1WithSpectralDensity(spectralDensityLiSz1 ...
        ,spectralDensityLiSz2,DD);
    meanR1WithLipariSzabo(atomNumber) = ...
        mean(r1WithLipariSzabo(1:atomNumber));
    
    r1WithPerturbationTheory(atomNumber) = ...
        calculateR1WithSpectralDensity(spectralDensity1 ...
        ,spectralDensity2,DD);
    meanR1WithPerturbationTheory(atomNumber) = mean( ...
        r1WithPerturbationTheory(1:atomNumber));
    
    if calculateSchroedingerEquation
        r1WithSchroedingerEquation(atomNumber) = ...
            calculateR1WithSchroedingerEquation(theta,phi ...
            ,nearestNeighbourDistancesPow3,deltaT,path2ConstantsFile);
        meanR1WithSchroedingerEquation(atomNumber) = ...
            mean(r1WithSchroedingerEquation(1:atomNumber));
    end
    
    if mod(atomNumber,20) == 0
        if calculateSchroedingerEquation
            save(path2Save ,'r1WithPerturbationTheory' ...
                ,'r1WithSchroedingerEquation' ...
                ,'meanR1WithPerturbationTheory' ...
                ,'meanR1WithSchroedingerEquation' ...
                ,'r1WithLipariSzabo' ...
                ,'meanR1WithLipariSzabo' ...
                ,'meanPositions','deltaT' ...
                ,'timeSteps','lags','numberOfHs','B0')
        else
            save(path2Save ,'r1WithPerturbationTheory' ...
                ,'meanR1WithPerturbationTheory' ...
                ,'r1WithLipariSzabo' ...
                ,'meanR1WithLipariSzabo' ...
                ,'meanPositions','deltaT' ...
                ,'timeSteps','lags','numberOfHs','B0')
        end
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

