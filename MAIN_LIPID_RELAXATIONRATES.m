%% load configuration
addpath('library')
configuration = readConfigurationFile('config.conf');

loaded = configuration.dataLoaded;
runOnServer = configuration.runOnServer; 

%% Load Trajectories
if not(loaded)
    
    clc
    clearvars -except configuration runOnServer
    close all
    
    loadReducedDataSet = configuration.loadReducedDataSet;
    
    LIPID = configuration.lipid;
    
    path2Project = configuration.path2Project;
    path2Data = configuration.path2Data;
    path2Results = configuration.path2Results;
    
    switch LIPID
        case{'DOPS'}
            fileName = configuration.fileNameDOPS;
            moleculeSize = configuration.moleculeSizeDOPS;
            headSize = configuration.headSizeDOPS;
        case{'PLPC'}
            fileName = configuration.fileNamePLPC;
            moleculeSize = configuration.moleculeSizePLPC;
            headSize = configuration.headSizePLPC;
        case{'PSM'}
            fileName = configuration.fileNamePSM;
            moleculeSize = configuration.moleculeSizePSM;
            headSize = configuration.headSizePSM;
        otherwise
            warining(['The lipid "' LIPD '" does not exist'])
    end
    
    path2Save = [path2Results fileName '_resultsCorrelationFunctions'];
    
    path2File = [path2Data fileName];
    if loadReducedDataSet
         path2File = [path2File '_reduced.mat'];
    else
         path2File = [path2File '.mat' ];
    end
     
    Data = load(path2File); 

    if loadReducedDataSet
        Data = Data(1).DataRed;
    end
     
    LipidH = Data.(fileName);

    trajX = squeeze(LipidH(:,1,:));
    trajY = squeeze(LipidH(:,2,:));
    trajZ = squeeze(LipidH(:,3,:));
end
%% Set Up
clearvars -except   deltaT ps trajX trajY trajZ path2Save moleculeSize ...
    LIPID headSize configuration runOnServer 
calculateSchroedingerEquation = configuration.calculateSchroedingerEquation;
path2ConstantsFile = configuration.path2ConstantsFile;
ps = 1e-12;
deltaT = configuration.deltaT*ps;

%% Define constants
constants = readConstantsFile(path2ConstantsFile);
hbar = constants.hbar;                              % [Js]
gammaRad = constants.gammaRad;                      % [rad/Ts]
B0 = constants.B0;                                  % [T]
mu0 = constants.mu0;                                % [N/A^2] Vacuum permeability
Nm = constants.Nm;                                  % [m] Nanometer (Traj. format)
DD = 3/4*(mu0/(4*pi)*hbar*gammaRad^2 )^2/(Nm^6);    % J/rad
omega0 = gammaRad*B0;                               % [rad/s]: Larmor (anglular) frequency

%% Define simulation parameters
[numberOfHs,timeSteps] = size(trajX)
lags = round(configuration.fractionForLags*timeSteps);                      
nearestNeighbours = configuration.nearestNeighbours;     
numberOfMolecules = round(numberOfHs/moleculeSize);


 %% Start simulation

r1WithPerturbationTheory = zeros(1,numberOfHs);
r1WithSchroedingerEquation = zeros(1,numberOfHs);

meanR1WithPerturbationTheory = zeros(1,numberOfHs);
meanR1WithSchroedingerEquation = zeros(1,numberOfHs);

meanPositions = zeros(numberOfHs,3);
molecules=1:numberOfHs;
%molecules=molecules(randperm(length(molecules)));

if ~runOnServer
    figure('Name','Relaxation Rates')
end

for atomNumber=1:numberOfHs
    tic
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
 
    % Szabo Fit
    
    
    
    [w1,w2] = calculateSpectralDensities(correlationFunction1, ...
        correlationFunction2,omega0,deltaT,lags);
    
    r1WithPerturbationTheory(atomNumber) = ...
        calculateR1WithSpectralDensity(w1(lags:end),w2(lags:end),DD);
    meanR1WithPerturbationTheory(atomNumber) = mean( ...
        r1WithPerturbationTheory(1:atomNumber));
    
    if calculateSchroedingerEquation
        r1WithSchroedingerEquation(atomNumber) = ...
            calculateR1WithSchroedingerEquation(theta,phi ...
            ,nearestNeighbourDistancesPow3,deltaT,path2ConstantsFile);
        meanR1WithSchroedingerEquation(atomNumber) = ...
            mean(r1WithSchroedingerEquation(1:atomNumber));
    end
    
    fileId = fopen('outputLog.txt','a');
    fprintf(fileId, ['Molecule: %f, Mean R1: \n' ...
        ,'   PerturbationTheory: %f \n' ...
        ,'   Schroedinger Equation: %f \n'],[atomNumber, ...
        r1WithPerturbationTheory(atomNumber) ...
        ,r1WithSchroedingerEquation(atomNumber)]);
    fclose(fileId);
    
    if mod(atomNumber,20)==0
        if calculateSchroedingerEquation
            save(path2Save ,'r1WithPerturbationTheory' ...
                ,'r1WithSchroedingerEquation' ...
                ,'meanR1WithPerturbationTheory' ...
                ,'meanR1WithSchroedingerEquation' ...
                ,'meanPositions','deltaT' ...
                ,'timeSteps','lags','numberOfHs','B0')
        else
            save(path2Save ,'r1WithPerturbationTheory' ...
                ,'meanR1WithPerturbationTheory' ...
                ,'meanPositions','deltaT' ...
                ,'timeSteps','lags','numberOfHs','B0')
        end
    end
    
    if ~runOnServer
        plot(r1WithPerturbationTheory(1:atomNumber),'b','LineWidth',1.5)
        hold on
        plot(meanR1WithPerturbationTheory(1:atomNumber),'--g' ... 
            ,'LineWidth',1.5)
        if calculateSchroedingerEquation
            plot(r1WithSchroedingerEquation(1:atomNumber),'r' ...
                ,'LineWidth',1.5)
            plot(meanR1WithSchroedingerEquation(1:atomNumber) ...
                ,'--y','LineWidth',1.5)
            legend('Spectral density','Mean Spectral Density' ...
                ,'Schroedinger Equation','Mean Schroedinger Equ.');
        else
            legend('Spectral density','Mean Spectral Density');
        end
        title('Relaxation Rate R1')
        xlabel('Epoches')
        ylabel('R1 [Hz]')
        grid on
        drawnow
        hold off
    end
    
end



    
 







