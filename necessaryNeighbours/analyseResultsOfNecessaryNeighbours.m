clc
clear all  %#ok<CLALL>
close all

results = load("C:\Users\maxoe\Google Drive\Promotion\Results\LIPIDS\PLPC\necessaryNearestNeighbours\20220425_Results_relevantNearestNeighbours_PLPClipid.mat");

r1 = results.r1WithPerturbationTheory;
whichLipid = results.whichLipid;
saving = 1;
savingPath = initializeSystemForSavingPlots("necessaryNearestNeighbours" ... 
    ,whichLipid);
atomCounter = results.atomCounter;

averageR1 = squeeze(mean(r1(:,:,1:atomCounter,:),3));
effectiveR1 = squeeze(mean(averageR1,2));
overallR1 = squeeze(mean(effectiveR1,1));

medianR1 = squeeze(median(r1(:,:,1:atomCounter,:),3));
effectiveMedianR1 = squeeze(mean(medianR1,2));
overallMedianR1 = squeeze(mean(effectiveMedianR1,1));

nearestNeighbourCases = results.nearestNeighbourCases;

rateShiftMean = effectiveR1 - effectiveR1(1,:);
rateShiftMedian = effectiveMedianR1 - effectiveMedianR1(1,:);

orientations = rad2deg(results.orientationAngles);

initializeFigure();
plot(nearestNeighbourCases,overallR1,'--')
plot(nearestNeighbourCases,overallMedianR1,'-.')
legend('Mean', 'Median','Location','East')
xlabel('Nearest Neighbours')
ylabel('Overall R$_1$ [Hz]')
title(sprintf('Overall nearest neighbour-dependent R$_1$ (%s)' ...
    ,whichLipid));

if saving
    savingName = sprintf('%s_overallR1NNDependent_%s' ...
        ,results.startDateOfSimulation,whichLipid);
    print(gcf,[savingPath savingName],'-dpng','-r300');
end

initializeFigure();
legendEntries = {};
for orientationNr = 1:length(orientations)
    plot(nearestNeighbourCases,effectiveR1(orientationNr,:),'*-')
    legendEntries{end+1} = sprintf('Mean, $\\theta$: %.2f ' ...
        ,orientations(orientationNr)); %#ok<SAGROW>
    plot(nearestNeighbourCases,effectiveMedianR1(orientationNr,:),'*-')
    legendEntries{end+1} = sprintf('Median, $\\theta$: %.2f' ...
        ,orientations(orientationNr)); %#ok<SAGROW>
end 
legend(legendEntries,'Location','East')
xlabel('Nearest neighbours')
ylabel('Relaxation rate [Hz]')
title(sprintf('Effective nearest neighbour-dependent R$_1$ (%s)' ...
    ,whichLipid));

if saving
    savingName = sprintf('%s_effectiveR1NNAndOrientationDependent_%s' ...
        ,results.startDateOfSimulation,whichLipid);
    print(gcf,[savingPath savingName],'-dpng','-r300');
end

initializeFigure();
legendEntries = {};
for orientationNr = 1:length(orientations)
    plot(nearestNeighbourCases,rateShiftMean(orientationNr,:),'*-')
    legendEntries{end+1} = sprintf('Mean, $\\theta$: %.2f ' ...
        ,orientations(orientationNr)); %#ok<SAGROW>
    plot(nearestNeighbourCases,rateShiftMedian(orientationNr,:),'*-')
    legendEntries{end+1} = sprintf('Median, $\\theta$: %.2f' ...
        ,orientations(orientationNr)); %#ok<SAGROW>
end 
xlabel('Nearest neighbours')
ylabel('Relaxation rate shift [Hz]')
title(sprintf('Nearest neighours and orientation-dependent R$_1$ shift (%s)' ...
    ,whichLipid));
legend(legendEntries,'Location','East')

if saving
    savingName = sprintf('%s_r1ShiftNNAndOrientationDependent_%s' ...
        ,results.startDateOfSimulation,whichLipid);
    print(gcf,[savingPath savingName],'-dpng','-r300');
end



