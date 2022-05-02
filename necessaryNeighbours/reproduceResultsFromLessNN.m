clc
clear all %#ok<CLALL>
close all

results = load("C:\Users\maxoe\Google Drive\Promotion\Results\LIPIDS\PLPC\necessaryNearestNeighbours\20220425_Results_relevantNearestNeighbours_PLPClipid.mat");
creationDate = results.startDateOfSimulation;
whichLipid = results.whichLipid;
savingPath = initializeSystemForSavingPlots('necessaryNearestNeighbours' ...
    ,whichLipid);
saving = 1;


r1 = results.r1WithPerturbationTheory;
orientationAngles = rad2deg(results.orientationAngles);
orientationsCount = length(orientationAngles);
positionAngles = rad2deg(results.positionAngles);
positionsCount = length(positionAngles);
calculatedAtoms = results.atomCounter;
nNCases = length(results.nearestNeighbourCases);

averagedR1 = squeeze(mean(r1,3));

differences = zeros(orientationsCount,positionsCount,nNCases);
sumDifferences = zeros(orientationsCount,positionsCount,nNCases);
summedRates = zeros(1,nNCases);
sumEstimatedR1s = zeros(orientationsCount,positionsCount,nNCases);
estimatedR1s = zeros(orientationsCount,positionsCount,nNCases);

for referenceOrientationNr = 1:orientationsCount
    for referencePositionNr = 1:positionsCount        
        referenceR1 = squeeze(averagedR1(referenceOrientationNr ...
            ,referencePositionNr,:))';
        highestNNR1 = referenceR1(1);
        rates = highestNNR1./referenceR1;
        summedRates = summedRates + rates;
        for orientationNr = 1:orientationsCount
            for positionNr = 1:positionsCount
                otherR1s = squeeze(averagedR1( ...
                    orientationNr,positionNr,:))';
                estimatedR1s(orientationNr,positionNr,:) ...
                    = otherR1s .* rates;
                differences(orientationNr,positionNr,:) = ...
                    otherR1s(1) ...
                - estimatedR1s(orientationNr,positionNr,:);
            end
        end
        averageEstimatedR1s(referenceOrientationNr ...
            ,referencePositionNr,:) = squeeze(mean(mean( ...
            estimatedR1s,2),1)); %#ok<SAGROW>
        sumDifferences(referenceOrientationNr ...
            ,referencePositionNr,:) = squeeze(sum(sum(abs( ...
            differences),2),1))';
    end
end

averageDifferences = sumDifferences/(orientationsCount*positionsCount-1);
averagedRates = summedRates/(orientationsCount*positionsCount);
averageEstimatedR1 = squeeze(mean(mean(averageEstimatedR1s,2),1));

sumDifferences = squeeze(mean(mean(abs(differences),2),1))';
nearestNeighbours = results.nearestNeighbourCases;

reproductionResults(1,:) = nearestNeighbours;
reproductionResults(2,:) = averagedRates;
reproductionResults(3,:) = squeeze(mean(mean(averageDifferences,2),1));
reproductionResults(4,:) = averageEstimatedR1;

initializeFigure();
legendEntries = {};
for referenceOrientationNr = 1:orientationsCount
    for referencePositionNr = 1:positionsCount
        plot(nearestNeighbours,squeeze(averageDifferences( ...
            referenceOrientationNr,referencePositionNr,:)),'*-')
        orientationAngle = orientationAngles(referenceOrientationNr);
        positionAngle = positionAngles(referencePositionNr);
        legendEntries{end+1} = sprintf("$\\theta$: %.1f $\\varphi$: %.1f" ...
            ,orientationAngle,positionAngle); %#ok<SAGROW>
    end
end
legend(legendEntries)
multiLineTitle = {sprintf('Reproduce R$_1$ with less nearest') ...
    , sprintf(' neighbours (Data: %s %s)',whichLipid,creationDate)};
title(multiLineTitle)
xlabel('Number of nearest neighbours to repoduce R$_1$')
ylabel('Averaged difference to original data')

if saving
    savingName = sprintf('%s_%s_%s.png',results.startDateOfSimulation ...
    ,'ReproduceRelaxationRates',whichLipid);
    print(gcf,[savingPath savingName],'-dpng','-r300');
    save(sprintf('%s%s.txt',savingPath,savingName),'reproductionResults' ...
        ,'-ASCII','-append');
end


