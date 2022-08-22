clc
clear all
close all

colorArray = [252, 3, 3;252, 194, 3;82, 252, 3;3, 252, 240;3, 252, 240; ...
    173, 3, 252;252, 3, 181;3, 186, 252;231, 252, 3;3, 152, 252]/252;
PSM = "C:\Users\maxoe\Google Drive\Promotion\Results\LIPIDS\PSM\orientationDependency\20220530_Results_orientationDependency_PSMlipid.mat";
DOPS = "C:\Users\maxoe\Google Drive\Promotion\Results\LIPIDS\DOPS\orientationDependency\20220307_Results_orientationDependency_DOPSlipid.mat";
PLPC = "C:\Users\maxoe\Google Drive\Promotion\Results\LIPIDS\PLPC\orientationDependency\20220502_Results_orientationDependency_PLPClipid.mat";
results = load(PLPC);
baseConfiguration = readConfigurationFile('../baseConfiguration.txt');
addpath(genpath(baseConfiguration.path2LibraryOnLocalMachine));

saving = 1;
if saving
    savingPath = initializeSystemForSavingPlots('reconstructionR1' ...
    ,results.whichLipid);
end

orientationAngles = rad2deg(results.orientationAngles);
positionAngles = rad2deg(results.positionAngles);

correlationFunctions1W0 = results.correlationFunction1W0Saver;
correlationFunctions2W0 = results.correlationFunction2W0Saver;
corrFuncLength = size(correlationFunctions1W0,3);

offsetReduction = 1;

if offsetReduction
    
    offsetReductionRegion = [round(0.5*corrFuncLength) ...
        :round(0.8*corrFuncLength)];
    initializeFigure('legend',false);
    title(sprintf('Unchanged correlation function (%s)', ...
        results.whichLipid));
    for orientationNr = 1:length(orientationAngles)
        for positionNr = 1:length(positionAngles)
            plot(squeeze(correlationFunctions1W0(orientationNr, ...
                positionNr,:)));
            plot(squeeze(correlationFunctions2W0(orientationNr, ...
                positionNr,:)));
        end
    end
    if saving
        savingName = sprintf('originalCorrFunctions');
        print(gcf,[savingPath savingName],'-dpng','-r300');
    end

    % subtract the offset
    for orientationNr = 1:length(orientationAngles)
        for positionNr = 1:length(positionAngles)
            correlationFunctions1W0(orientationNr,positionNr,:) = ...
                correlationFunctions1W0(orientationNr,positionNr,:) ...
                - mean(correlationFunctions1W0(orientationNr,positionNr ...
                ,offsetReductionRegion));
            correlationFunctions2W0(orientationNr,positionNr,:) = ...
                correlationFunctions2W0(orientationNr,positionNr,:) ...
                - mean(correlationFunctions2W0(orientationNr,positionNr ...
                ,offsetReductionRegion));
        end
    end
end

% cut correlation function to avoid the increase of the correlation
% function caused by FFT
cutOffFraction = 0.7;
cutOffRegion = [1:round(corrFuncLength*cutOffFraction)-1];
for orientationNr = 1:length(orientationAngles)
    for positionNr = 1:length(positionAngles)
        cuttedCorrelationFunctions1W0(orientationNr,positionNr,:) = ...
            correlationFunctions1W0(orientationNr,positionNr,cutOffRegion);
        cuttedCorrelationFunctions2W0(orientationNr,positionNr,:) = ...
            correlationFunctions2W0(orientationNr,positionNr,cutOffRegion);
    end
end

% plot correlation functions after cut and offset substraction
initializeFigure('legend',false);
title(sprintf('Offset substracted and cutted correlation functions (%s)', ...
    results.whichLipid));
for orientationNr = 1:length(orientationAngles)
    for positionNr = 1:length(positionAngles)
        plot(squeeze(cuttedCorrelationFunctions1W0(orientationNr ...
            ,positionNr,:)));
        plot(squeeze(cuttedCorrelationFunctions2W0(orientationNr ...
            ,positionNr,:)));
    end
end
if saving
    savingName = sprintf('offsetReductedAndCuttedCorrFunctions');
    print(gcf,[savingPath savingName],'-dpng','-r300');
end

for orientationNr = 1:length(orientationAngles)
        for positionNr = 1:length(positionAngles)
            correlationFunctions1W0(orientationNr,positionNr,:) = ...
                correlationFunctions1W0(orientationNr,positionNr,:) ...
                - mean(correlationFunctions1W0(orientationNr,positionNr ...
                ,offsetReductionRegion));
            correlationFunctions2W0(orientationNr,positionNr,:) = ...
                correlationFunctions2W0(orientationNr,positionNr,:) ...
                - mean(correlationFunctions2W0(orientationNr,positionNr ...
                ,offsetReductionRegion));
        end
end

timeSkips = [1 2 3 5 10];
for timeNr = 1:length(timeSkips)
    timeSkip = timeSkips(timeNr);
    deltaT = results.samplingFrequency*timeSkip;
    for orientationNr = 1:length(orientationAngles)
        for positionNr = 1:length(positionAngles)
            correlationFunction1W0 = squeeze(cuttedCorrelationFunctions1W0( ...
                orientationNr,positionNr,1:timeSkip:end))';
            correlationFunction2W0 = squeeze(cuttedCorrelationFunctions2W0( ...
                orientationNr,positionNr,1:timeSkip:end))';
            [spectralDensity1W0,spectralDensity2W0] = ...
                calculateSpectralDensities(correlationFunction1W0 ...
                ,correlationFunction2W0,results.omega0 ...
                ,deltaT,length(correlationFunction1W0));
            r1WithPerturbationTheory(timeNr,orientationNr,positionNr) = ...
                calculateR1WithSpectralDensity(spectralDensity1W0 ...
                ,spectralDensity2W0,results.dipolDipolConstant) ...
                *results.configuration.scalingRate;  %#ok<SAGROW>
            
        end
    end
end


r1Sim = squeeze(mean(results.scaledR1Rates(:,:,1:results.atomCounter),3));

lineStyles = [":" "--" "-." "-"]; 
initializeFigure();
legendEntries = {};
for positionNr = 1:length(positionAngles)
    plot(orientationAngles,r1WithPerturbationTheory(1, ...
        :,positionNr),lineStyles(1),'color',colorArray(positionNr,:));
    legendEntries{end+1} = sprintf('Reconstructed, %.2f', ...
        positionAngles(positionNr)); %#ok<SAGROW>
    plot(orientationAngles,r1Sim(:,positionNr) ...
        ,lineStyles(2),'color',colorArray(positionNr,:));
    legendEntries{end+1} = sprintf('Simulated, %.2f', ...
        positionAngles(positionNr));  %#ok<SAGROW>

end
legend(legendEntries);
xlabel(sprintf('Orientation Angle $\\theta$'));
ylabel(sprintf('Relaxation rate [Hz]'));
 title(sprintf('Reconstructed vs. simulated R1 (%s)', ...
        results.whichLipid));
if saving
    savingName = sprintf('reconstructedVsSimulatedR1');
    print(gcf,[savingPath savingName],'-dpng','-r300');
end

%% print results
fprintf('%s:\n',results.whichLipid);
fprintf('  averaged R1 Simulation:\n    deltaT = %.2d: %.4f \n', ...
    results.samplingFrequency,mean(mean(r1Sim)));
fprintf('  averaged R1 Reconstruction:\n')
for timeNr = 1:length(timeSkips)
    deltaT = results.samplingFrequency*timeSkips(timeNr);
    fprintf('    deltaT = %.2d: %.4f\n', ...
        deltaT,mean(mean(r1WithPerturbationTheory(timeNr,:,:),3),2));
end

% fprintf(['%s: \n    averaged relaxation rate simulation: %.4f \n' ...
%     '    averaged relaxation rate reconstruction: %.4f\n'], ...
%     results.whichLipid,mean(mean(r1Sim)), ...
%     mean(mean(r1WithPerturbationTheory,3),2));




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



