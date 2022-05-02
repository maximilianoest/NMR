clc
clear all
close all

results = load("C:\Users\maxoe\Google Drive\Promotion\Results\locationDependentR1\20220426_Results_locationDependentR1_DOPSlipid4.mat");
baseConfiguration = readConfigurationFile('../baseConfiguration.txt');
addpath(genpath(baseConfiguration.path2LibraryOnLocalMachine));

saving = 1;
savingPath = sprintf(['C:\\Users\\maxoe\\Google Drive\\Promotion' ...
    '\\Zwischenergebnisse\\%s_simulationLength_%s\\'] ...
    ,datestr(date,'yyyymmdd'),results.whichLipid);
if ~exist(savingPath,'dir')
    mkdir(savingPath);
end

orientationAngles = rad2deg(results.orientationAngles);
samplingFrequency = results.samplingFrequency;
lineStyle = ["-","--","-.",":"];
lineStyle = [lineStyle lineStyle lineStyle lineStyle lineStyle lineStyle];


fractionNames = fieldnames(results.correlationFunction0W0Saver);
maxSize = size(results.correlationFunction0W0Saver.(fractionNames{1}),3);
timeAxis = [0:maxSize-1] * samplingFrequency;
for orientationNr = 1:length(orientationAngles)
    currentFigure = initializeFigure();
    legendEntries = {};
    orientationAngle = orientationAngles(orientationNr);
    for simTimeNr = 1:2:length(fractionNames)
        simTimeFraction = results.simulationTimeFractions(simTimeNr);
        fractionName = fractionNames{simTimeNr};
        corrFunc = results.correlationFunction0W0Saver.(fractionName);
        effCorrFunc = squeeze(mean(corrFunc(orientationNr,:,:),2));
        normEffCorrFunc = effCorrFunc'/effCorrFunc(1);
        normEffCorrFunc(end+1:maxSize+1) = NaN;
        plot(timeAxis,normEffCorrFunc(1:end-1),lineStyle(simTimeNr))
        legendEntries{end+1} = sprintf('SimTimeFraction: %.2f' ...
            ,simTimeFraction); %#ok<SAGROW>
    end
    xlabel('Correlation time [s]')
    axis([0 inf -inf 0.2])
    title(sprintf(['Normalized correlation 0 $\\omega_0$ functions ' ...
        'for $\\theta$: %.2f $^\\circ$'],orientationAngle));
    legend(legendEntries)
    
    initializeFigure();
    legendEntries = {};
    orientationAngle = orientationAngles(orientationNr);
    for simTimeNr = 1:2:length(fractionNames)
        simTimeFraction = results.simulationTimeFractions(simTimeNr);
        fractionName = fractionNames{simTimeNr};
        corrFunc = results.correlationFunction1W0Saver.(fractionName);
        effCorrFunc = squeeze(mean(corrFunc(orientationNr,:,:),2));
        normEffCorrFunc = effCorrFunc'/effCorrFunc(1);
        normEffCorrFunc(end+1:maxSize+1) = NaN;
        plot(timeAxis,normEffCorrFunc(1:end-1),lineStyle(simTimeNr))
        legendEntries{end+1} = sprintf('SimTimeFraction: %.2f' ...
            ,simTimeFraction); %#ok<SAGROW>
    end
    xlabel('Correlation time [s]')
    axis([0 inf -inf 0.2])
    title(sprintf(['Normalized correlation 1 $\\omega_0$ functions ' ...
        'for $\\theta$: %.2f $^\\circ$'],orientationAngle));
    legend(legendEntries)
    
    initializeFigure();
    legendEntries = {};
    orientationAngle = orientationAngles(orientationNr);
    for simTimeNr = 1:2:length(fractionNames)
        simTimeFraction = results.simulationTimeFractions(simTimeNr);
        fractionName = fractionNames{simTimeNr};
        corrFunc = results.correlationFunction2W0Saver.(fractionName);
        effCorrFunc = squeeze(mean(corrFunc(orientationNr,:,:),2));
        normEffCorrFunc = effCorrFunc'/effCorrFunc(1);
        normEffCorrFunc(end+1:maxSize+1) = NaN;
        plot(timeAxis,normEffCorrFunc(1:end-1),lineStyle(simTimeNr))
        legendEntries{end+1} = sprintf('SimTimeFraction: %.2f' ...
            ,simTimeFraction); %#ok<SAGROW>
    end
    xlabel('Correlation time [s]')
    axis([0 inf -inf 0.2])
    title(sprintf(['Normalized correlation 2 $\\omega_0$ functions ' ...
        'for $\\theta$: %.2f $^\\circ$'],orientationAngle));
    legend(legendEntries)
    
end
% atomCounter = results.atomCounter;
% effCorrelationFunction0W0 = squeeze(mean(results.sumCorrelationFunction0W0Saver,2))/atomCounter;
% effCorrelationFunction1W0 = squeeze(mean(results.sumCorrelationFunction1W0Saver,2))/atomCounter;
% effCorrelationFunction2W0 = squeeze(mean(results.sumCorrelationFunction2W0Saver,2))/atomCounter;
% 
% effCorrelationFunction0W0 = effCorrelationFunction0W0 ...
%     ./ effCorrelationFunction0W0(:,1);
% effCorrelationFunction1W0 = effCorrelationFunction1W0 ...
%     ./ effCorrelationFunction1W0(:,1);
% effCorrelationFunction2W0 = effCorrelationFunction2W0 ...
%     ./ effCorrelationFunction2W0(:,1);
% 
% orientationAngles = results.orientationAngles;
% samplingFrequency = results.samplingFrequency;
% timeAxis = [0:length(effCorrelationFunction0W0)-1] * samplingFrequency;
% 
% 
% for orientationNr = 1:length(orientationAngles)
%     figure;
%     hold on
%     plot(timeAxis,abs(effCorrelationFunction0W0(orientationNr,:)) ...
%         ,'LineWidth',1.5)
%     plot(timeAxis,abs(effCorrelationFunction1W0(orientationNr,:)) ...
%         ,'LineWidth',1.5)
%     plot(timeAxis,abs(effCorrelationFunction2W0(orientationNr,:)) ...
%         ,'LineWidth',1.5)
%     hold off
%     grid minor
%     legend('0 \omega_0','1 \omega_0','2 \omega_0')
%     title(['Normalized N correlation functions for orientation ' ...
%         num2str(rad2deg(orientationAngles(orientationNr))) '°']) 
%     xlabel('Correlation time [s]')
%     axis([0 inf 0 0.2])
%     
% end





% %% load data
% clc
% clear all
% close all
% 
% 
% configuration = readConfigurationFile('config.conf');
% addpath(genpath(configuration.path2Library));
% compartment = configuration.compartment;
% 
% data = loadResultsFromR1Simulation(configuration);
% 
% atomCount = data.atomCounter;
% orientationAngles = rad2deg(data.orientationAngles);
% orientationsCount = length(orientationAngles);
% positionAngles = rad2deg(data.positionAngles);
% positionsCount = length(positionAngles);
% deltaT = data.deltaT*data.shiftForCorrelationFunction;
% 
% B0 = data.B0;
% B0WithoutComma = manipulateB0ValueForSavingPath(B0);
% 
% startDateOfSimulation = data.startDateOfSimulation;
% path2SaveFigures = [configuration.path2Results startDateOfSimulation ...
%     '_CorrelationFunctionsOf' compartment 'At' ...
%     num2str(B0WithoutComma) 'T.fig'];
% 
% %% analyze data
% [correlationFunctions1W0,correlationFunctions2W0] = ...
%     getCorrelationFunctionsFromSimulationData(data);
% 
% averageCorrelationFunctions1W0 = squeeze(mean( ...
%     correlationFunctions1W0,3));
% averageCorrelationFunctions2W0 = squeeze(mean( ...
%     correlationFunctions2W0,3));
% 
% effectiveCorrelationFunctions1W0 = squeeze(mean( ...
%     averageCorrelationFunctions1W0,2));
% effectiveCorrelationFunctions2W0 = squeeze(mean( ...
%     averageCorrelationFunctions2W0,2));
% 
% %% check for difference in phi
% differenceForPosition1W0 = zeros(orientationsCount,positionsCount-1, ...
%     size(effectiveCorrelationFunctions1W0,2));
% differenceForPosition2W0 = zeros(orientationsCount,positionsCount-1, ...
%     size(effectiveCorrelationFunctions1W0,2));
% 
% for orientationNumber = 1:orientationsCount
%     for positionNumber = 1:positionsCount-1
%         differenceForPosition1W0(orientationNumber,positionNumber,:) = ...
%             abs(real(squeeze(averageCorrelationFunctions1W0( ...
%             orientationNumber,1,:)-averageCorrelationFunctions1W0( ...
%             orientationNumber,positionNumber+1,:))'));
%         differenceForPosition2W0(orientationNumber,positionNumber,:) = ...
%             abs(real(squeeze(averageCorrelationFunctions2W0( ...
%             orientationNumber,1,:)-averageCorrelationFunctions2W0( ...
%             orientationNumber,positionNumber+1,:))'));
%     end
% end
% 
% %% plotting
% tauAxis = 0:deltaT:(size(correlationFunctions1W0,4)-1)*deltaT;
% tauMin = 0;
% tauMax = 5e-7;
% valueMin = 0;
% valueMax = 0.5;
% fontSize = 15;
% 
% legendEntries = {};
% figs(1) = figure('DefaultAxesFontSize',fontSize);
% hold on
% for orientationNumber = 1:orientationsCount
%     plot(tauAxis,abs(real(effectiveCorrelationFunctions1W0( ...
%         orientationNumber,:)/effectiveCorrelationFunctions1W0( ...
%         orientationNumber,1))),'LineWidth',1.5)
%     legendEntries{orientationNumber} = ['\theta: ' num2str( ...
%         orientationAngles(orientationNumber))]; %#ok<SAGROW>
% end
% axis([tauMin tauMax valueMin valueMax])
% hold off
% grid minor
% legend(legendEntries)
% title('Effective Correlation Functions for w0')
% xlabel('\tau')
% 
% legendEntries = {};
% figs(2) = figure('DefaultAxesFontSize',fontSize);
% hold on
% legendEntryCounter = 1;
% for orientationNumber = 1:orientationsCount
%     for positionNumber = 1:positionsCount
%         plot(tauAxis,squeeze(abs(real(averageCorrelationFunctions1W0( ...
%             orientationNumber,positionNumber,:) ...
%             /averageCorrelationFunctions1W0( ...
%             orientationNumber,positionNumber,1)))),'LineWidth',1.5)
%         legendEntries{legendEntryCounter} = ['\theta: ' ...
%             num2str(orientationAngles(orientationNumber)) ', \phi: ' ...
%             num2str(positionAngles(positionNumber))]; %#ok<SAGROW>
%         legendEntryCounter = legendEntryCounter+1;
%     end
% end
% axis([tauMin tauMax valueMin valueMax])
% hold off
% grid minor
% legend(legendEntries)
% title('Average Correlation Functions for w0')
% xlabel('\tau')
% 
% legendEntries = {};
% figs(3) = figure('DefaultAxesFontSize',fontSize);
% hold on
% for orientationNumber = 1:orientationsCount
%     plot(tauAxis,abs(real(effectiveCorrelationFunctions2W0( ...
%         orientationNumber,:)/effectiveCorrelationFunctions2W0( ...
%         orientationNumber,1))),'LineWidth',1.5)
%     legendEntries{orientationNumber} = ['\theta: ' num2str( ...
%         orientationAngles(orientationNumber))]; %#ok<SAGROW>
%     
% end
% axis([tauMin tauMax valueMin valueMax])
% hold off
% grid minor
% legend(legendEntries)
% title('Effective Correlation Functions for 2w0')
% xlabel('\tau')
% 
% legendEntries = {};
% figs(4) = figure('DefaultAxesFontSize',fontSize);
% hold on
% legendEntryCounter = 1;
% for orientationNumber = 1:orientationsCount
%     for positionNumber = 1:positionsCount
%         plot(tauAxis,squeeze(abs(real(averageCorrelationFunctions2W0( ...
%             orientationNumber,positionNumber,:) ...
%             /averageCorrelationFunctions2W0( ...
%             orientationNumber,positionNumber,1)))),'LineWidth',1.5)
%         legendEntries{legendEntryCounter} = ['\theta: ' ...
%             num2str(orientationAngles(orientationNumber)) ', \phi: ' ...
%             num2str(positionAngles(positionNumber))]; %#ok<SAGROW>
%         legendEntryCounter = legendEntryCounter+1;
%     end
%     
% end
% axis([tauMin tauMax valueMin valueMax])
% hold off
% grid minor
% legend(legendEntries)
% title('Average Correlation Functions for w0')
% xlabel('\tau')
% ylabel('G(\tau)')
% 
% legendEntries = {};
% fig(5) = figure('DefaultAxesFontSize',fontSize);
% hold on
% legendEntryCounter = 1;
% for orientationNumber = 1:orientationsCount
%     orientationAngle = orientationAngles(orientationNumber);
%     for positionNumber = 1:positionsCount-1
%         positionAngle = positionAngles(positionNumber+1);
%         
%         legendEntries{legendEntryCounter} = ['\omega_0' ...
%             '\theta: ' num2str(orientationAngle) ', ' ...
%             '\phi: ' num2str(positionAngle)];  %#ok<SAGROW>
%         legendEntryCounter = legendEntryCounter+1;
%         plot(tauAxis,squeeze(differenceForPosition1W0(orientationNumber ...
%             ,positionNumber,:)),'LineWidth',1.5)
%         
%         legendEntries{legendEntryCounter} = ['2\omega_0' ...
%             '\theta: ' num2str(orientationAngle) ', ' ...
%             '\phi: ' num2str(positionAngle)]; %#ok<SAGROW>
%         legendEntryCounter = legendEntryCounter+1;
%         plot(tauAxis,squeeze(differenceForPosition2W0(orientationNumber ...
%             ,positionNumber,:)),'LineWidth',1.5)
%     end
% end
% hold off
% ylabel(['Difference to Position at: ' num2str(positionAngles(1)) 'ï¿½'])
% xlabel('\tau')
% legend(legendEntries)
% grid minor
% 
% savefig(figs,path2SaveFigures);
% 
% 
% 
% 
