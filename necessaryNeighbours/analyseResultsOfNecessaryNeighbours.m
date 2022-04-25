clc
clear all  %#ok<CLALL>
close all

% configuration = readConfigurationFile('config.txt');
% if configuration.runOnServer
%     path2Results = configuration.path2ResultsOnServer;
%     addpath(genpath(configuration.path2LibraryOnServer));
% else
%     path2Results = configuration.path2ResultsOnLocalMachine;
%     addpath(genpath(configuration.path2LibraryOnLocalMachine));
% end
results = load("C:\Users\maxoe\Google Drive\Promotion\Results\LIPIDS\PLPC\necessaryNearestNeighbours\20220419_Results_relevantNearestNeighbours_PLPClipid.mat");

r1 = results.r1WithPerturbationTheory;
saving = 1;
savingPath = sprintf(['C:\\Users\\maxoe\\Google Drive\\Promotion' ...
    '\\Zwischenergebnisse\\%s_numberOfAtoms_%s\\'] ...
    ,datestr(date,'yyyymmdd'),results.whichLipid);
if ~exist(savingPath,'dir')
    mkdir(savingPath);
end
atomCounter = results.atomCounter;

averageR1 = squeeze(mean(r1(:,:,1:atomCounter,:),3));
effectiveR1 = squeeze(mean(averageR1,2));
overallR1 = squeeze(mean(effectiveR1,1));

medianR1 = squeeze(median(r1(:,:,1:atomCounter,:),3));
effectiveMedianR1 = squeeze(mean(medianR1,2));
overallMedianR1 = squeeze(mean(effectiveMedianR1,1));

whichLipid = results.whichLipid;
nearestNeighbourCases = results.nearestNeighbourCases;

rateShiftMean = effectiveR1 - effectiveR1(1,:);
rateShiftMedian = effectiveMedianR1 - effectiveMedianR1(1,:);

orientations = rad2deg(results.orientationAngles);

figurePosAndSize = [50 50 900 600];
figure('Position',figurePosAndSize);
set(gcf,'DefaultLineLineWidth',1.5)
set(gca,'FontSize',16)
hold on
plot(nearestNeighbourCases,overallR1,'--','LineWidth', 1.5)
plot(nearestNeighbourCases,overallMedianR1,'-.','LineWidth',1.5)
hold off
legend('Mean', 'Median','Location','East')
xlabel('Nearest Neighbours')
ylabel('Overall R1 [Hz]')
title(sprintf('Overall nearest neighbour-dependent R_1 (%s)' ...
    ,whichLipid));
grid minor

if saving
    savingName = sprintf('overallR1NNDependent_%s',whichLipid);
    print(gcf,[savingPath savingName],'-dpng','-r300');
end

figurePosAndSize = [50 50 900 600];
figure('Position',figurePosAndSize);
set(gcf,'DefaultLineLineWidth',1.5)
set(gca,'FontSize',16)
hold on
legendEntries = {};
for orientationNr = 1:length(orientations)
    plot(nearestNeighbourCases,effectiveR1(orientationNr,:),'*-' ...
        ,'LineWidth', 1.5)
    legendEntries{end+1} = sprintf('Mean, \\theta: %.2f ' ...
        ,orientations(orientationNr)); %#ok<SAGROW>
    plot(nearestNeighbourCases,effectiveMedianR1(orientationNr,:) ...
        ,'*-','LineWidth', 1.5)
    legendEntries{end+1} = sprintf('Median, \\theta: %.2f' ...
        ,orientations(orientationNr)); %#ok<SAGROW>
end 
hold off
legend(legendEntries,'Location','East')
xlabel('Nearest neighbours')
ylabel('Relaxation rate [Hz]')
title(sprintf('Effective nearest neighbour-dependent R_1 (%s)' ...
    ,whichLipid));
grid minor

if saving
    savingName = sprintf('effectiveR1NNAndOrientationDependent_%s',whichLipid);
    print(gcf,[savingPath savingName],'-dpng','-r300');
end

figurePosAndSize = [50 50 900 600];
figure('Position',figurePosAndSize);
set(gcf,'DefaultLineLineWidth',1.5)
set(gca,'FontSize',16)
hold on
legendEntries = {};
for orientationNr = 1:length(orientations)
    plot(nearestNeighbourCases,rateShiftMean(orientationNr,:) ...
        ,'*-','LineWidth', 1.5)
    legendEntries{end+1} = sprintf('Mean, \\theta: %.2f ' ...
        ,orientations(orientationNr)); %#ok<SAGROW>
    plot(nearestNeighbourCases,rateShiftMedian(orientationNr,:) ...
        ,'*-','LineWidth', 1.5)
    legendEntries{end+1} = sprintf('Median, \\theta: %.2f' ...
        ,orientations(orientationNr)); %#ok<SAGROW>
end 
hold off
xlabel('Nearest neighbours')
ylabel('Relaxation rate shift [Hz]')
title(sprintf('Nearest neighours and orientation-dependent R_1 (%s)' ...
    ,whichLipid));
grid minor
legend(legendEntries,'Location','East')

if saving
    savingName = sprintf('r1ShiftNNAndOrientationDependent_%s',whichLipid);
    print(gcf,[savingPath savingName],'-dpng','-r300');
end



