clc
clear all  %#ok<CLALL>
close all

PSM = "C:\Users\maxoe\Google Drive\Promotion\Results\LIPIDS\PSM\necessaryNearestNeighbours\20220517_Results_relevantNearestNeighbours_PSMlipid.mat";
DOPS = "FREE";
PLPC = "C:\Users\maxoe\Google Drive\Promotion\Results\LIPIDS\PLPC\necessaryNearestNeighbours\20220425_Results_relevantNearestNeighbours_PLPClipid.mat";
results = load(PLPC);

r1 = results.r1WithPerturbationTheory;
whichLipid = results.whichLipid;
configuration = readConfigurationFile("config.txt");
baseConfiguration =  readConfigurationFile(configuration ...
    .path2BaseConfigurationOnLocalMachine);
addpath(genpath(baseConfiguration.path2LibraryOnLocalMachine))
creationDate = results.startDateOfSimulation;
nearestNeighbours = results.nearestNeighbourCases;
constants = readConstantsFile("../../constants.txt");
omega0 = constants.gyromagneticRatioOfHydrogenAtom ...
    *3;%0.276;
saving = 0;
plotCorrelationFunctions = 1;
if saving
    savingPath = initializeSystemForSavingPlots("necessaryNearestNeighbours" ...
        ,whichLipid);
end
atomCounter = results.atomCounter;

orientationAngles = rad2deg(results.orientationAngles);
orientationsCount = length(orientationAngles);
positionAngles = rad2deg(results.positionAngles);
positionsCount = length(positionAngles);

correlationFunctions1W0 = results.correlationFunction1W0Saver;
correlationFunctions2W0 = results.correlationFunction2W0Saver;
corrFuncLength = size(correlationFunctions1W0,4);
nearestNeighbourCases = results.nearestNeighbourCases;
nNCases = length(nearestNeighbourCases);

offsetSuppressionFractions = [0.5 0.8];
cutOffFraction = 0.7;

% subtract offset of correlation function resulting from time constant
% Hamiltonian part.
offsetSuppressionRegion = round(offsetSuppressionFractions(1) ...
    *corrFuncLength):round(offsetSuppressionFractions(2)*corrFuncLength);
for nnCaseNr = 1:length(nearestNeighbourCases)
    for orientationNr = 1:length(orientationAngles)
        for positionNr = 1:length(positionAngles)
            correlationFunctions1W0(nnCaseNr,orientationNr,positionNr,:) ...
                = correlationFunctions1W0(nnCaseNr,orientationNr,positionNr,:) ...
                - mean(correlationFunctions1W0(nnCaseNr,orientationNr,positionNr ...
                ,offsetSuppressionRegion));
            correlationFunctions2W0(nnCaseNr,orientationNr,positionNr,:) ...
                = correlationFunctions2W0(nnCaseNr,orientationNr,positionNr,:) ...
                - mean(correlationFunctions2W0(nnCaseNr,orientationNr,positionNr ...
                ,offsetSuppressionRegion));
        end
    end
end

% cut correlation function to avoid the increase of the correlation
% function at higher tau caused by FFT.
cutOffRegion = [1:round(corrFuncLength*cutOffFraction)-1];
for nnCaseNr = 1:length(nearestNeighbourCases)
    for orientationNr = 1:length(orientationAngles)
        for positionNr = 1:length(positionAngles)
            cuttedCorrelationFunctions1W0(nnCaseNr,orientationNr, ...
                positionNr,:) = correlationFunctions1W0(nnCaseNr, ...
                orientationNr,positionNr,cutOffRegion);
            cuttedCorrelationFunctions2W0(nnCaseNr,orientationNr, ...
                positionNr,:) = correlationFunctions2W0(nnCaseNr, ...
                orientationNr,positionNr,cutOffRegion);
        end
    end
end

% plot correlation functions after cut and offset substraction
if plotCorrelationFunctions
    for nnCaseNr = 1:2:length(nearestNeighbourCases)
        initializeFigure('legend',false);
        title(sprintf(['Offset substracted and cutted correlation functions (NN:%i %s)'], ...
            nearestNeighbourCases(nnCaseNr),results.whichLipid));
        for orientationNr = 1:length(orientationAngles)
            for positionNr = 1:length(positionAngles)
                plot(squeeze(cuttedCorrelationFunctions1W0(nnCaseNr, ...
                    orientationNr,positionNr,:)));
                plot(squeeze(cuttedCorrelationFunctions2W0(nnCaseNr, ...
                    orientationNr,positionNr,:)));
            end
        end
        if saving
            savingName = sprintf('%s_offsetSupprAndCuttedCorrFunc_NN%i', ...
                results.startDateOfSimulation, ...
                nearestNeighbourCases(nnCaseNr));
            print(gcf,[savingPath savingName],'-dpng','-r300');
        end
    end
end


deltaT = results.samplingFrequency;
for nnCaseNr = 1:length(nearestNeighbourCases)
    for orientationNr = 1:length(orientationAngles)
        for positionNr = 1:length(positionAngles)
            correlationFunction1W0 = squeeze(cuttedCorrelationFunctions1W0( ...
                nnCaseNr,orientationNr,positionNr,:))';
            correlationFunction2W0 = squeeze(cuttedCorrelationFunctions2W0( ...
                nnCaseNr,orientationNr,positionNr,:))';
            [spectralDensity1W0,spectralDensity2W0] = ...
                calculateSpectralDensities(correlationFunction1W0 ...
                ,correlationFunction2W0,omega0 ...
                ,deltaT,length(correlationFunction1W0));
            avgNNR1(orientationNr,positionNr,nnCaseNr) = ...
                calculateR1WithSpectralDensity(spectralDensity1W0 ...
                ,spectralDensity2W0,results.dipolDipolConstant); %#ok<SAGROW>
        end
    end
end

% calculate scaling factor for 300 atoms and validate it on 100 atoms ->
% not possible because the correlation function is averaged over the 400
% atoms.

% 1. determine scaling rate based on calculated R1 in reconstruction
effNNR1 = squeeze(mean(avgNNR1,2));
overallNNR1 = squeeze(mean(effNNR1,1));

initializeFigure('legend', false);
plot(nearestNeighbourCases,overallNNR1,'*-')
xlabel('Nearest Neighbours')
ylabel('Overall R$_1$ [Hz]')
title(sprintf('Overall R$_1$ (%s)' ...
    ,whichLipid));

if saving
    savingName = sprintf('%s_overallR1NNDependent_%s' ...
        ,results.startDateOfSimulation,whichLipid);
    print(gcf,[savingPath savingName],'-dpng','-r300');
end

initializeFigure();
legendEntries = {};
for orientationNr = 1:length(orientationAngles)
    plot(nearestNeighbourCases,effNNR1(orientationNr,:),'*-')
    legendEntries{end+1} = sprintf('Mean, $\\theta$: %.2f ' ...
        ,orientationAngles(orientationNr)); %#ok<SAGROW>
end
legend(legendEntries,'Location','East')
xlabel('Nearest neighbours')
ylabel('R$_1$[Hz]')
title(sprintf('Effective R$_1$ (%s)' ...
    ,whichLipid));

if saving
    savingName = sprintf('%s_effectiveR1NNDependent_%s' ...
        ,results.startDateOfSimulation,whichLipid);
    print(gcf,[savingPath savingName],'-dpng','-r300');
end

initializeFigure();
legendEntries = {};
for nnCaseNr = 1:length(nearestNeighbourCases)
   plot(orientationAngles,effNNR1(:,nnCaseNr));
   legendEntries{end+1} = sprintf('%d NN',nearestNeighbourCases(nnCaseNr));
end
legend(legendEntries,'Location','northwest');
xlabel('Orientation angle $\theta$');
ylabel('R$_1$ [Hz]');
title(sprintf('Effective R$_1$ (%s)',whichLipid));

if saving
    savingName = sprintf('%s_effectiveR1ThetaDependent_%s' ...
        ,results.startDateOfSimulation,whichLipid);
    print(gcf,[savingPath savingName],'-dpng','-r300');
end

initializeFigure();
legendEntries = {};
thetaDependentR1Shift = effNNR1 - effNNR1(1,:);

for orientationNr = 1:length(orientationAngles)
    plot(nearestNeighbourCases,thetaDependentR1Shift( ...
        orientationNr,:),'*-');
    legendEntries{end+1} = sprintf('$\\theta$: %.2f ' ...
        ,orientationAngles(orientationNr)); %#ok<SAGROW>
end
xlabel('Nearest neighbours')
ylabel('R$_1$ shift [Hz]')
title(sprintf('$\\theta$-dependent R$_1$ shift (%s)' ...
    ,whichLipid));
legend(legendEntries,'Location','East')

if saving
    savingName = sprintf('%s_effectiveR1ShiftNNDependent_%s' ...
        ,results.startDateOfSimulation,whichLipid);
    print(gcf,[savingPath savingName],'-dpng','-r300');
end

%% REFERENCE SCALING: predict and plot results from less NN 
predictedR1FromLessNN = zeros(orientationsCount,positionsCount, ...
    orientationsCount,positionsCount,nNCases);
scalingRatesAngleDependent = zeros(orientationsCount,positionsCount,nNCases);

% 1. predict R1 (theta and phi dependent)
for refOriNr = 1:orientationsCount
    for refPosNr = 1:positionsCount
        referenceR1 = squeeze(avgNNR1(refOriNr ...
            ,refPosNr,:))';
        highestNNR1 = referenceR1(1);
        scalingRates = highestNNR1./referenceR1;
        scalingRatesAngleDependent(refOriNr, ...
            refPosNr,:) = scalingRates;
        for orientationNr = 1:orientationsCount
            for positionNr = 1:positionsCount
                allNNR1 = squeeze(avgNNR1( ...
                    orientationNr,positionNr,:))';
                predictedR1FromLessNN(refOriNr, ...
                    refPosNr,orientationNr,positionNr,:) ...
                    = allNNR1 .* scalingRates;
            end
        end
    end
end

% 2. plot theta-dependent R1 values
thetaDepPredictedR1 = squeeze(mean(predictedR1FromLessNN,4));
initializeFigure();
legendEntries = {};

for refOriNr = 1:size(thetaDepPredictedR1,1)
    for refPosNr = 1:size(thetaDepPredictedR1,2)
        for oriNr = 1:size(thetaDepPredictedR1,3)
            p = plot(nearestNeighbourCases, ...
                squeeze(thetaDepPredictedR1(refOriNr,refPosNr,oriNr,:)));
            if oriNr == 1
                oldColor = p.Color;
            else
                p.Color = oldColor;
            end
            legendEntries{end+1} = sprintf( ...
                '$\\theta$: %4.2f (Ref.: $\\theta$: %4.2f, $\\varphi$: %4.2f)', ...
                orientationAngles(oriNr),orientationAngles(refOriNr),positionAngles(refPosNr));
        end
    end
end

legend(legendEntries,'location','east');
xlabel('Nearest neighbours')
ylabel('R$_1$ [Hz]')
title(sprintf('$\\theta$-dependent predicted R$_1$ (%s)' ...
    ,whichLipid));
legend(legendEntries,'Location','East')

if saving
    savingName = sprintf('%s_predictedR1_thetaNNDependent_%s' ...
        ,results.startDateOfSimulation,whichLipid);
    print(gcf,[savingPath savingName],'-dpng','-r300');
end

% 3. plot theta-dependent R1 shift
initializeFigure();
legendEntries={};

for refOriNr = 1:orientationsCount
    for refPosNr = 1:positionsCount
        predictedR1 = squeeze(thetaDepPredictedR1(refOriNr,refPosNr,:,:,:));
        scaledThetatDependentR1Shift = predictedR1 - predictedR1(1,:);
        for oriNr = 1:orientationsCount
            p = plot(nearestNeighbourCases, ...
                    scaledThetatDependentR1Shift(oriNr,:));
            if oriNr == 1
                oldColor = p.Color;
            else
                p.Color = oldColor;
            end
            legendEntries{end+1} = sprintf('$\\theta$: %4.2f (Ref.: $\\theta$: %4.2f, $\\varphi$: %4.2f)',orientationAngles(oriNr),orientationAngles(refOriNr),positionAngles(refPosNr));
        end
    end
end
legend(legendEntries,'location','east')
xlabel('Nearest neighbours')
ylabel('R$_1$ shift [Hz]')
title(sprintf('$\\theta$-dependent predicted R$_1$ shift (%s)' ...
    ,whichLipid));
legend(legendEntries,'Location','East')

if saving
    savingName = sprintf('%s_predictedR1Shift_thetaNNDependent_%s' ...
        ,results.startDateOfSimulation,whichLipid);
    print(gcf,[savingPath savingName],'-dpng','-r300');
end


% 4. plot predicited overall R1
overallPredictedR1 = squeeze(mean(thetaDepPredictedR1,3));

initializeFigure();
legendEntries = {};
for refOriNr = 1:orientationsCount
    for refPosNr = 1:positionsCount
        plot(nearestNeighbourCases,squeeze(overallPredictedR1(refOriNr,refPosNr,:)),'*-');
        legendEntries{end+1} = sprintf('Ref.: $\\theta$: %.2f, $\\varphi$: %.2f',orientationAngles(refOriNr),positionAngles(refPosNr));
    end
end

xlabel('Nearest neighbours')
ylabel('Predicted R$_1$ [Hz]')
title(sprintf('NN-dependent predicted R$_1$  (%s)',whichLipid));
legend(legendEntries,'Location','East')

if saving
    savingName = sprintf('%s_predictedR1_NNDependent_%s' ...
        ,results.startDateOfSimulation,whichLipid);
    print(gcf,[savingPath savingName],'-dpng','-r300');
end

% 5. plot difference in overall R1 to highest NN case
overallDifference = overallPredictedR1 - overallPredictedR1(:,:,1);

initializeFigure();
legendEntries = {};
for refOriNr = 1:orientationsCount
    for refPosNr = 1:positionsCount
        plot(nearestNeighbourCases,squeeze(overallDifference(refOriNr,refPosNr,:)),'*-');
        legendEntries{end+1} = sprintf('Ref.: $\\theta$: %.2f, $\\varphi$: %.2f',orientationAngles(refOriNr),positionAngles(refPosNr));
    end
end
legend(legendEntries,'location','east')

xlabel('Nearest neighbours')
ylabel('Difference [Hz]')
title(sprintf('Difference in predicted overall R$_1$ to higher NN (%s)',whichLipid));

if saving
    savingName = sprintf('%s_predictedR1DiffrenceToHigh_NNDependent_%s' ...
        ,results.startDateOfSimulation,whichLipid);
    print(gcf,[savingPath savingName],'-dpng','-r300');
end

%% determine theta- and phi-dependent R1 with avg scaling rate

%% determine theta-dependent R1 predicted with average scaling rate
avgScalingRates = mean(mean(scalingRatesAngleDependent,2),1);
predictedR1WithAvgScalRate = avgNNR1.*avgScalingRates;
effPredictedR1WithAvgScalRate = squeeze(mean(predictedR1WithAvgScalRate,2));

initializeFigure();
legendEntries = {};
for oriNr = 1:size(thetaDepPredictedR1,3)
    p = plot(nearestNeighbourCases, ...
        squeeze(effPredictedR1WithAvgScalRate(oriNr,:)));
    legendEntries{end+1} = sprintf('$\\theta$: %4.2f',orientationAngles(oriNr));
end
legend(legendEntries);
xlabel('Nearest neighbours')
ylabel('R$_1$ [Hz]')
title(sprintf('$\\theta$-dependent predicted R$_1$ with averaged scaling rates (%s)' ...
    ,whichLipid));

if saving
    savingName = sprintf('%s_effectiveAvgScalingPredictedR1NNDependent_%s' ...
        ,results.startDateOfSimulation,whichLipid);
    print(gcf,[savingPath savingName],'-dpng','-r300');
end

%%
averageDifferences = sumDifferences/(orientationsCount*positionsCount-1);
averagedScalingRates = summedScalingRates/(orientationsCount*positionsCount);
averagepredicteddR1 = squeeze(mean(mean(overallPredictedR1s,2),1));

initializeFigure();
legendEntries = {};
for refOriNr = 1:orientationsCount
    for refPosNr = 1:positionsCount
        plot(nearestNeighbours,squeeze(averageDifferences( ...
            refOriNr,refPosNr,:)),'*-')
        orientationAngle = orientationAngles(refOriNr);
        positionAngle = positionAngles(refPosNr);
        legendEntries{end+1} = sprintf("Ref. $\\theta$: %.1f $\\varphi$: %.1f" ...
            ,orientationAngle,positionAngle); %#ok<SAGROW>
    end
end
legend(legendEntries)
title(sprintf('Predicted reconstructed R$_1$ with scaling rate (Data: %s %s)', ...
    whichLipid,creationDate));
xlabel('Number of nearest neighbours to predict R$_1$')
ylabel('Averaged difference to higher NN')





reconstructedScaledResults(1,:) = nearestNeighbours;
reconstructedScaledResults(2,:) = averagedScalingRates;
reconstructedScaledResults(3,:) = squeeze(mean(mean(averageDifferences,2),1));
reconstructedScaledResults(4,:) = averagepredicteddR1;
reconstructedScaledResults(5,:) = overallNNR1;

if saving
    savingName = sprintf('%s_%s_%s.png',results.startDateOfSimulation ...
        ,'ReproduceRelaxationRates',whichLipid);
    print(gcf,[savingPath savingName],'-dpng','-r300');
    save(sprintf('%s%s.txt',savingPath,savingName),'scaledResults' ...
        ,'-ASCII','-append');
end


% 2. determine scaling rate based on calcuated R1 in simulation
%%
% avgSimNNR1 = squeeze(mean(r1(:,:,1:atomCounter,:),3));
% effSimNNR1 = squeeze(mean(avgSimNNR1,2));
% overallSimNNR1 = squeeze(mean(effSimNNR1,1));
% 
% medSimNNR1 = squeeze(median(r1(:,:,1:atomCounter,:),3));
% effMedSimNNR1 = squeeze(mean(medSimNNR1,2));
% overallMedSimNNR1 = squeeze(mean(effMedSimNNR1,1));
% 
% rateDifferenceMean = effSimNNR1 - effSimNNR1(1,:);
% rateDifferenceMedian = effMedSimNNR1 - effMedSimNNR1(1,:);
% 
% initializeFigure();
% plot(nearestNeighbourCases,overallSimNNR1,'-*')
% % plot(nearestNeighbourCases,overallMedSimNNR1,'-.')
% legend('Mean', 'Median','Location','East')
% xlabel('Nearest Neighbours')
% ylabel('Overall R$_1$ [Hz]')
% title(sprintf('Overall nearest neighbour-dependent R$_1$ (%s)' ...
%     ,whichLipid));
% 
% if saving
%     savingName = sprintf('%s_overallR1NNDependent_%s' ...
%         ,results.startDateOfSimulation,whichLipid);
%     print(gcf,[savingPath savingName],'-dpng','-r300');
% end
% 
% %%
% initializeFigure();
% legendEntries = {};
% for orientationNr = 1:length(orientationAngles)
%     plot(nearestNeighbourCases,effSimNNR1(orientationNr,:),'*-')
%     legendEntries{end+1} = sprintf('Mean, $\\theta$: %.2f ' ...
%         ,orientationAngles(orientationNr)); %#ok<SAGROW>
%     plot(nearestNeighbourCases,effMedSimNNR1(orientationNr,:),'*-')
%     legendEntries{end+1} = sprintf('Median, $\\theta$: %.2f' ...
%         ,orientationAngles(orientationNr)); %#ok<SAGROW>
% end
% legend(legendEntries,'Location','East')
% xlabel('Nearest neighbours')
% ylabel('Relaxation rate [Hz]')
% title(sprintf('Effective nearest neighbour-dependent R$_1$ (%s)' ...
%     ,whichLipid));
% 
% if saving
%     savingName = sprintf('%s_effectiveR1NNAndOrientationDependent_%s' ...
%         ,results.startDateOfSimulation,whichLipid);
%     print(gcf,[savingPath savingName],'-dpng','-r300');
% end
% 
% initializeFigure();
% legendEntries = {};
% for orientationNr = 1:length(orientationAngles)
%     plot(nearestNeighbourCases,rateDifferenceMean(orientationNr,:),'*-')
%     legendEntries{end+1} = sprintf('Mean, $\\theta$: %.2f ' ...
%         ,orientationAngles(orientationNr)); %#ok<SAGROW>
%     plot(nearestNeighbourCases,rateDifferenceMedian(orientationNr,:),'*-')
%     legendEntries{end+1} = sprintf('Median, $\\theta$: %.2f' ...
%         ,orientationAngles(orientationNr)); %#ok<SAGROW>
% end
% xlabel('Nearest neighbours')
% ylabel('Relaxation rate difference to highest NN [Hz]')
% title(sprintf('Nearest neighours and orientation-dependent R$_1$ shift (%s)' ...
%     ,whichLipid));
% legend(legendEntries,'Location','East')
% 
% if saving
%     savingName = sprintf('%s_differenceR1NNAndOrientationDependent_%s' ...
%         ,results.startDateOfSimulation,whichLipid);
%     print(gcf,[savingPath savingName],'-dpng','-r300');
% end
% 
% differences = zeros(orientationsCount,positionsCount,nNCases);
% sumDifferences = zeros(orientationsCount,positionsCount,nNCases);
% summedScalingRates = zeros(1,nNCases);
% sumEstimatedR1s = zeros(orientationsCount,positionsCount,nNCases);
% predictedR1s = zeros(orientationsCount,positionsCount,nNCases);
% 
% for referenceOrientationNr = 1:orientationsCount
%     for referencePositionNr = 1:positionsCount
%         referenceR1 = squeeze(avgSimNNR1(referenceOrientationNr ...
%             ,referencePositionNr,:))';
%         highestNNR1 = referenceR1(1);
%         scalingRates = highestNNR1./referenceR1;
%         summedScalingRates = summedScalingRates + scalingRates;
%         for orientationNr = 1:orientationsCount
%             for positionNr = 1:positionsCount
%                 r1sForAllNN = squeeze(avgSimNNR1( ...
%                     orientationNr,positionNr,:))';
%                 predictedR1s(orientationNr,positionNr,:) ...
%                     = r1sForAllNN .* scalingRates;
%                 differences(orientationNr,positionNr,:) = ...
%                     r1sForAllNN(1) ...
%                     - predictedR1s(orientationNr,positionNr,:);
%             end
%         end
%         overallPredictedR1s(referenceOrientationNr ...
%             ,referencePositionNr,:) = squeeze(mean(mean( ...
%             predictedR1s,2),1)); %#ok<SAGROW>
%         sumDifferences(referenceOrientationNr ...
%             ,referencePositionNr,:) = squeeze(sum(sum(abs( ...
%             differences),2),1))';
%     end
% end
% 
% averageDifferences = sumDifferences/(orientationsCount*positionsCount-1);
% averagedScalingRates = summedScalingRates/(orientationsCount*positionsCount);
% averagepredicteddR1 = squeeze(mean(mean(overallPredictedR1s,2),1));
% 
% sumDifferences = squeeze(mean(mean(abs(differences),2),1))';
% nearestNeighbours = results.nearestNeighbourCases;
% 
% scaledResults(1,:) = nearestNeighbours;
% scaledResults(2,:) = averagedScalingRates;
% scaledResults(3,:) = squeeze(mean(mean(averageDifferences,2),1));
% scaledResults(4,:) = averagepredicteddR1;
% 
% initializeFigure();
% legendEntries = {};
% for referenceOrientationNr = 1:orientationsCount
%     for referencePositionNr = 1:positionsCount
%         plot(nearestNeighbours,squeeze(averageDifferences( ...
%             referenceOrientationNr,referencePositionNr,:)),'*-')
%         orientationAngle = orientationAngles(referenceOrientationNr);
%         positionAngle = positionAngles(referencePositionNr);
%         legendEntries{end+1} = sprintf("Ref. $\\theta$: %.1f $\\varphi$: %.1f" ...
%             ,orientationAngle,positionAngle); %#ok<SAGROW>
%     end
% end
% legend(legendEntries)
% title(sprintf('Predicted R$_1$ with scaling rate (Data: %s %s)', ...
%     whichLipid,creationDate));
% xlabel('Number of nearest neighbours to repoduce R$_1$')
% ylabel('Averaged difference to higher NN')
% 
% if saving
%     savingName = sprintf('%s_%s_%s.png',results.startDateOfSimulation ...
%         ,'ReproduceRelaxationRates',whichLipid);
%     print(gcf,[savingPath savingName],'-dpng','-r300');
%     save(sprintf('%s%s.txt',savingPath,savingName),'scaledResults' ...
%         ,'-ASCII','-append');
% end



