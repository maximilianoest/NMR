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
ps = 1e-12;
deltaT = configuration.deltaT*ps;


%% Define constants
constants = readConstantsFile('constants.txt');
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
offset = configuration.offsetFromHead;
snippetSize = configuration.snippetSize;  
numberOfMolecules = round(numberOfHs/moleculeSize);

%% Analyse Trajectories and Isolate molecule snippet

atomIndices = getIndicesOfSnippet(moleculeSize,headSize,offset ...
    ,snippetSize, numberOfMolecules);

if ~runOnServer
    
    figure('Name','All Molecules')
    width = 4;

    hold on
    for moleculeNumber = 0:numberOfMolecules
        firstIndex = moleculeNumber*moleculeSize+1;
        if (moleculeNumber+1)*moleculeSize > numberOfHs
            break
        else
            secondIndex = (moleculeNumber+1)*moleculeSize;
        end
        plot3(trajX(firstIndex:secondIndex,end) ...
            ,trajY(firstIndex:secondIndex,end) ...
            ,trajZ(firstIndex:secondIndex,end),'-x','Linewidth',width);
    end
    hold off
    xlabel('X')
    ylabel('Y')
    zlabel('Z')
    legend('1', '2', '3')
    grid on

    cuttedRegionsX = NaN(1,numberOfHs);
    cuttedRegionsY = NaN(1,numberOfHs);
    cuttedRegionsZ = NaN(1,numberOfHs);

    for moleculeNr = 1:numberOfMolecules
        for atomNumber = 1:10
            cuttedRegionsX(atomIndices(moleculeNr,atomNumber)) ...
                = trajX(atomIndices(moleculeNr,atomNumber),end);
            cuttedRegionsY(atomIndices(moleculeNr,atomNumber)) ...
                = trajY(atomIndices(moleculeNr,atomNumber),end);
            cuttedRegionsZ(atomIndices(moleculeNr,atomNumber)) ...
                = trajZ(atomIndices(moleculeNr,atomNumber),end);
        end
    end

    figure('Name','Isolated Parts of Molecules')
    width = 4;
    hold on
    for moleculeNumber = 0:numberOfMolecules
        firstIndex = moleculeNumber*moleculeSize+1;
        if (moleculeNumber+1)*moleculeSize > numberOfHs
            break
        else
            secondIndex = (moleculeNumber+1)*moleculeSize;
        end
        plot3(trajX(firstIndex:secondIndex,end) ...
            ,trajY(firstIndex:secondIndex,end) ...
            ,trajZ(firstIndex:secondIndex,end) ...
            ,'-x','Linewidth',width);
        plot3(cuttedRegionsX(firstIndex:secondIndex) ...
            ,cuttedRegionsY(firstIndex:secondIndex) ...
            ,cuttedRegionsZ(firstIndex:secondIndex) ...
            ,'-og','LineWidth',width);
    end
    hold off
    xlabel('X')
    ylabel('Y')
    zlabel('Z')
    grid on
end


%% calculate correlation function

correlation = zeros(numberOfMolecules,snippetSize,2,lags);

for moleculeNumber = 1:numberOfMolecules
    disp(sprintf('Molecule %i ',moleculeNumber))
    for atomNumber=1:snippetSize
        disp(sprintf('Atom %i ',atomNumber))
        indx = atomIndices(moleculeNumber,atomNumber);
        
        [relativeXPositions,relativeYPositions,relativeZPositions ...
            ,distances] = calculatePositionsAndDistances(trajX,trajY ...
            ,trajZ,indx);

        [nearestNeighboursX,nearestNeighboursY,nearestNeighboursZ ...
            ,nearestNeighbourDistancesPow3] = findNearestNeighbours( ...
            distances,nearestNeighbours+1,indx,relativeXPositions...
            ,relativeYPositions,relativeZPositions);  
        
        [F1,F2] = calculateSphericalHarmonics(nearestNeighboursX ...
            ,nearestNeighboursY,nearestNeighboursZ ...
            ,nearestNeighbourDistancesPow3);

        correlationFunction1 = calculateCrossCorrelationFunction(F1 ...
            ,conj(F1),lags-1);
        correlationFunction2 = calculateCrossCorrelationFunction(F2 ... 
            ,conj(F2),lags-1);
        correlation(moleculeNumber,atomNumber,1,:) ...
            = mean(correlationFunction1,1);
        correlation(moleculeNumber,atomNumber,2,:) ...
            = mean(correlationFunction2,1);
    end
    
end

%% plot correlation function
meanCorrelation = squeeze(mean(squeeze(mean(abs(real(correlation)),2)),1));
normCorrelation(1,:) = squeeze(meanCorrelation(1,:)) ...
    /(meanCorrelation(1,1));
normCorrelation(2,:) = squeeze(meanCorrelation(2,:)) ...
    /(meanCorrelation(2,1));

save(path2Save,'meanCorrelation','normCorrelation');
if ~runOnServer
    timeAxis = 0:deltaT:(lags-1)*deltaT;
    figure('Name','Correlation Functions') 
    plot(timeAxis,normCorrelation(1,:),'LineWidth',1.2)
    hold on
    plot(timeAxis,normCorrelation(2,:),'LineWidth',1.2)
    hold off
    grid on
    legend('w_0', '2 x w_0')
    xlabel('tau')
    title(['Correlation Function of ',LIPID])
    axis([0 timeAxis(end) 0 0.5])
end

 %% Start simulation
% 
% r1SpectralDensity = zeros(1,numberOfHs);
% r1SchroedingerEquation = zeros(1,numberOfHs);
% 
% meanR1SpectralDensity = zeros(1,numberOfHs);
% meanR1SchroedingerEquation = zeros(1,numberOfHs);
% 
% meanPositions = zeros(numberOfHs,3);
% 
% molecules=1:numberOfHs;
% %molecules=molecules(randperm(length(molecules)));
% 
% figure('Name','Relaxation Rates')
% 
% for atomNumber=1:numberOfHs
%     tic
%     
%     % positions of obeserved H in world coordinate system
%     xPosition = trajX(molecules(atomNumber),:);
%     yPosition = trajY(molecules(atomNumber),:);
%     zPosition = trajZ(molecules(atomNumber),:);
%     
% %     % mean position of oberserved H in time
% %     meanPositions(moleculeNr,:)=[mean(xPosition) ...
% %         ,mean(yPosition),mean(zPosition)];
%     
%     % positions of other particles relative to observed H
%     relativeXPositions = bsxfun(@minus,trajX,xPosition);
%     relativeYPositions = bsxfun(@minus,trajY,yPosition);
%     relativeZPositions = bsxfun(@minus,trajZ,zPosition);
%     distance = sqrt(relativeXPositions.^2+relativeYPositions.^2 ... 
%         +relativeZPositions.^2);
%     
%     % Find the n nearest neighbours 
%     inverseMeanDistance = 1./mean(distance,2);
%     [~,indx,~] = findNearestNeighbours(inverseMeanDistance...
%         ,nearestNeighbours+1);
%     id = false(size(inverseMeanDistance));
%     id(indx) = true;
%     id(molecules(atomNumber)) = false;
% 
%     nearestNeighboursX = relativeXPositions(id,:);
%     nearestNeighboursY = relativeYPositions(id,:);
%     nearestNeighboursZ = relativeZPositions(id,:);
%     nearestNeighbourDistancesPow3 = distance(id,:).^3;    
% 
%     % cartesian to spherical coordinates
%     hypotuseXY = hypot(nearestNeighboursX,nearestNeighboursY);
%     theta =pi/2-atan2(nearestNeighboursZ,hypotuseXY);
%     phi = atan2(nearestNeighboursY,nearestNeighboursX);
%         
%     SIN = sin(theta);
%     COS = cos(theta);
%          
%     % spherical harmonics
%     F1 = SIN.*COS.*exp(-1i*phi)./nearestNeighbourDistancesPow3;
%     F2 = SIN.^2.*exp(-2i*phi)./nearestNeighbourDistancesPow3;
%         
%     % correlation functions
%     correlationFunction1 = crosscorrelationOverMultipleDimensions(F1 ...
%         ,conj(F1),lags-1);
%     correlationFunction2 = crosscorrelationOverMultipleDimensions(F2 ... 
%         ,conj(F2),lags-1);
%     
%     % spectral densities
%     w1 = 2*(deltaT*cumsum(sum(correlationFunction1) ... 
%         .*exp(-1i*omega0*deltaT*(0:lags-1))))  ;
%     w2 = 2*(deltaT*cumsum(sum(correlationFunction2) ...
%         .*exp(-2i*omega0*deltaT*(0:lags-1))))  ;
%     
%     % calculate relaxation rate by using spectral density
%     r1SpectralDensity(atomNumber) = DD*3/2*(abs(real(mean(w1(lags:end)))) ... 
%         +abs(real(mean(w2(lags:end)))));
%     meanR1SpectralDensity(atomNumber) = mean(r1SpectralDensity( ...
%         1:atomNumber));
%     
%     % calculate relaxation rate by solving Schroedingers Equation
%     if calculateSchroedingerEquation
%         cur_R1SchroedingerEquation = getR1FromSchroedingerEquation(B0 ...
%             ,theta,phi,nearestNeighbourDistancesPow3,deltaT);
%         r1SchroedingerEquation(atomNumber) = ...
%             sum(cur_R1SchroedingerEquation);
%         meanR1SchroedingerEquation(atomNumber) = ...
%             mean(r1SchroedingerEquation(1:atomNumber));
%     end
%     
%     disp(sprintf('Molecule %i, Mean R1:',atomNumber))
%     disp(sprintf('  spectral density: %f',mean(r1SpectralDensity)))
%     if calculateSchroedingerEquation
%         disp(sprintf('  schroedinger equation: %f' ...
%             ,meanR1SchroedingerEquation))
%     end
%     
% %     SumCor1(moleculeNr,:)=sum(correlationFunction1);
% %     SumCor2(moleculeNr,:)=sum(correlationFunction2);
%     
%     plot(r1SpectralDensity(1:atomNumber),'b','LineWidth',1.5)
%     hold on
%     plot(meanR1SpectralDensity(1:atomNumber),'--g','LineWidth',1.5)
%     if calculateSchroedingerEquation
%         plot(r1SchroedingerEquation(1:atomNumber),'r','LineWidth',1.5)
%         plot(meanR1SchroedingerEquation(1:atomNumber),'--y','LineWidth',1.5)
%         legend('Spectral density','Mean Spectral Density' ...
%             ,'Schroedinger Equation','Mean Schroedinger Equ.');
%     else
%         legend('Spectral density','Mean Spectral Density');
%     end
%     title('Relaxation Rate R1')
%     xlabel('Epoches')
%     ylabel('R1 [Hz]')
%     grid on
%     drawnow
%     hold off
%     
%     if mod(atomNumber,20)==0
%         if calculateSchroedingerEquation
%          save(Path2Save ,'meanR1SpectralDensity' ...
%              ,'meanr1SchroedingerEquation','MeanPositions','DeltaT' ...
%              ,'timeSteps','lags','numberOfHs','B0')
%         else
%             save(Path2Save ,'meanR1SpectralDensity' ...
%              ,'MeanPositions','DeltaT','timeSteps','lags','numberOfHs' ...
%              ,'B0')
%         end
%     end
%     
% end



    
 







