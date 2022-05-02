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
positionAngles = rad2deg(results.positionAngles);

correlationFunctions1W0 = results.correlationFunction1W0Saver;
correlationFunctions2W0 = results.correlationFunction2W0Saver;

fractionNames = results.fractionNames;

for timeNr = 1:length(fractionNames)
    fractionName = fractionNames{timeNr};
    for orientationNr = 1:length(orientationAngles)
        for positionNr = 1:length(positionAngles)
            correlationFunction1W0 = squeeze(correlationFunctions1W0.( ...
                fractionName)(orientationNr,positionNr,:))';
            correlationFunction2W0 = squeeze(correlationFunctions2W0.( ...
                fractionName)(orientationNr,positionNr,:))';
            [spectralDensity1W0,spectralDensity2W0] = ...
                calculateSpectralDensities(correlationFunction1W0 ...
                ,correlationFunction2W0,results.omega0 ...
                ,results.samplingFrequency,length(correlationFunction1W0));
            r1WithPerturbationTheory(orientationNr,positionNr,timeNr) = ...
                calculateR1WithSpectralDensity(spectralDensity1W0 ...
                ,spectralDensity2W0,results.dipolDipolConstant);  %#ok<SAGROW>
            
        end
    end
end
r1Sim = squeeze(mean(results.r1WithPerturbationTheory,3));
lineStyles = [":" "--" "-." "-"]; 
for positionNr = 1:length(positionAngles)
    legendEntries = {};
    initializeFigure('legendFontSize',12);
    
    lineStyleCounter = 1;
    for timeNr = 1:2:length(fractionNames)
        plot(orientationAngles,r1WithPerturbationTheory( ...
            :,positionNr,timeNr),lineStyles(lineStyleCounter));
        legendEntries{end+1} = sprintf('Reconstruction, %.2f' ...
            ,results.simulationTimeFractions(timeNr)); %#ok<SAGROW>
        
        plot(orientationAngles,r1Sim(:,positionNr,timeNr) ...
            ,lineStyles(lineStyleCounter));
        legendEntries{end+1} = sprintf('Simulation, %.2f' ...
            ,results.simulationTimeFractions(timeNr));  %#ok<SAGROW>
        
        lineStyleCounter = lineStyleCounter+1;
    end
    title(sprintf('Position $\\varphi$: %.2f',positionAngles(positionNr)));
    xlabel(sprintf('Orientation Angle $\\theta$'));
    ylabel(sprintf('Relaxation rate (%s)',results.whichLipid));
    legend(legendEntries);
    if saving
        savingName = sprintf('simulationLengthDependentR1_phi%i_%s' ...
            ,positionAngles(positionNr),results.whichLipid);
        print(gcf,[savingPath savingName],'-dpng','-r300');
    end

end





% disp('-- Only summed up correlation function --')
% corrFunc1W0 = results.sumCorrelationFunction1W0Saver/results.atomCounter;
% corrFunc2W0 = results.sumCorrelationFunction2W0Saver/results.atomCounter;
% 
% for orientationNr = 1:length(orientationAngles)
%     for positionNr = 1:length(positionAngles)
%         [spectralDensity1W0, spectralDensity2W0] = ...
%             calculateSpectralDensities(squeeze(corrFunc1W0(orientationNr ...
%             ,positionNr,:))',squeeze(corrFunc2W0(orientationNr,positionNr,:))' ...
%             ,results.omega0,results.samplingFrequency,results.lags);
%         r1WithPerturbationTheory(orientationNr,positionNr) = ...
%             calculateR1WithSpectralDensity(spectralDensity1W0 ...
%             ,spectralDensity2W0,results.dipolDipolConstant);
%     end
% end
% r1WithPerturbationTheory
% r1 = mean(results.r1WithPerturbationTheory(:,:,:,1),3)

% corrFunc1W0 = results.correlationFunction1W0Saver;
% corrFunc2W0 = results.correlationFunction2W0Saver;
% 
% disp('-- Averaging already in simulation --')
% 
% for orientationNr = 1:length(orientationAngles)
%     for positionNr = 1:length(positionAngles)
%         [spectralDensity1W0, spectralDensity2W0] = ...
%             calculateSpectralDensities(squeeze(corrFunc1W0(orientationNr ...
%             ,positionNr,:))',squeeze(corrFunc2W0(orientationNr,positionNr,:))' ...
%             ,results.omega0,results.samplingFrequency,results.lags);
%         r1WithPerturbationTheory(orientationNr,positionNr) = ...
%             calculateR1WithSpectralDensity(spectralDensity1W0 ...
%             ,spectralDensity2W0,results.dipolDipolConstant);
%     end
% end
% r1WithPerturbationTheory
% r1 = mean(results.r1WithPerturbationTheory(:,:,:),3)


% corrFunc1W0 = results.allCorrelationFunction1W0Saver;
% corrFunc2W0 = results.allCorrelationFunction2W0Saver;
% 
% disp('-- First caclulate R1 then average R1 --')
% for atomCounter = 1:results.atomCounter
%     for orientationNr = 1:length(orientationAngles)
%         for positionNr = 1:length(positionAngles)
%             
%             [spectralDensity1W0, spectralDensity2W0] = ...
%                 calculateSpectralDensities(squeeze(corrFunc1W0(atomCounter,orientationNr ...
%                 ,positionNr,:))',squeeze(corrFunc2W0(atomCounter,orientationNr,positionNr,:))' ...
%                 ,results.omega0,results.samplingFrequency,results.lags);
%             r1WithPerturbationTheory(atomCounter,orientationNr,positionNr) = ...
%                 calculateR1WithSpectralDensity(spectralDensity1W0 ...
%                 ,spectralDensity2W0,results.dipolDipolConstant);
%         end
%     end
% end
% 
% squeeze(mean(r1WithPerturbationTheory,1))
% r1 = mean(results.r1WithPerturbationTheory(:,:,:,1),3)



