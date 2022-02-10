clc
clear all  %#ok<CLALL>
close all

configuration = readConfigurationFile('config.txt');
if configuration.runOnServer
    path2Results = configuration.path2ResultsOnServer;
    addpath(genpath(configuration.path2LibraryOnServer));
else
    path2Results = configuration.path2ResultsOnLocalMachine;
    addpath(genpath(configuration.path2LibraryOnLocalMachine));
end

load([path2Results configuration.fileName '.mat']);

relaxationRates = r1WithPerturbationTheory; 
atomCounter = atomCounter;

averageR1 = squeeze(mean(relaxationRates(:,:,1:atomCounter,:),3));
effectiveR1 = squeeze(mean(averageR1,2));
overallR1 = squeeze(mean(effectiveR1,1));

medianR1 = squeeze(median(relaxationRates(:,:,1:atomCounter,:),3));
effectiveMedianR1 = squeeze(mean(medianR1,2));
overallMedianR1 = squeeze(mean(effectiveMedianR1,1));

fileName = strsplit(fileName,'_');
fileName = fileName{1};

rateShiftMean = effectiveR1 - effectiveR1(1,:);
rateShiftMedian = effectiveMedianR1 - effectiveMedianR1(1,:);

orientations = rad2deg(orientationAngles);
material = strsplit(fileName,'_');
material = material{1};

figure(1)
hold on
plot(nearestNeighbourCases,overallR1,'--','LineWidth', 1.5)
plot(nearestNeighbourCases,overallMedianR1,'-.','LineWidth',1.5)
hold off
legend('Mean', 'Median')
xlabel('Nearest Neighbours')
ylabel('Overall R1 [Hz]')
title(['Necessity of number of nearest Neighbours (' material ')'])
grid minor
creationDate = datestr(now,'yyyymmdd');
filePath = sprintf('%s%s_%s_%s.png',path2Results,creationDate ...
    ,material,'EffRelaxationRateNNDependent');
saveas(gcf,filePath);

figure(2)
legendEntries = {};
hold on
for orientationCounter = 1:length(orientations)
    plot(nearestNeighbourCases,effectiveR1(orientationCounter,:),'*-' ...
        ,'LineWidth', 1.5)
    legendEntries{end+1} = ['Mean, Orientation ' ...
        num2str(orientations(orientationCounter))];
    plot(nearestNeighbourCases,effectiveMedianR1(orientationCounter,:) ...
        ,'*-','LineWidth', 1.5)
    legendEntries{end+1} = ['Median, Orientation ' ...
        num2str(orientations(orientationCounter))];
end 
hold off
legend(legendEntries)
xlabel('Nearest neighbours')
ylabel('Relaxation rate [Hz]')
title(['Nearest Neighbours for orientation dependency (' material ')'])
grid minor
creationDate = datestr(now,'yyyymmdd');
filePath = sprintf('%s%s_%s_%s.png',path2Results,creationDate ...
    ,material,'RelaxationRateNNDependent');
saveas(gcf,filePath);

figure(3)
legendEntries = {};
hold on
for orientationCounter = 1:length(orientations)
    plot(nearestNeighbourCases,rateShiftMean(orientationCounter,:) ...
        ,'*-','LineWidth', 1.5)
    legendEntries{end+1} = ['Mean, Orientation ' ...
        num2str(orientations(orientationCounter))];
    plot(nearestNeighbourCases,rateShiftMedian(orientationCounter,:) ...
        ,'*-','LineWidth', 1.5)
    legendEntries{end+1} = ['Median, Orientation ' ...
        num2str(orientations(orientationCounter))];
end 
hold off
xlabel('Nearest neighbours')
ylabel('Relaxation rate shift [Hz]')
title(['Rate shifts in dependence of nearest neighbours and' ...
    ' orientation (' material ')'])
grid minor
legend(legendEntries)
creationDate = datestr(now,'yyyymmdd');
filePath = sprintf('%s%s_%s_%s.png',path2Results,creationDate ...
    ,material,'orientationShiftNNDependent');
saveas(gcf,filePath);


