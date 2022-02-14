clc
clear all %#ok<CLALL>
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
creationDate = strsplit(configuration.fileName,'_');
creationDate = creationDate{1};

r1Results = r1WithPerturbationTheory;
[orientationsCount,positionsCount,calculatedAtoms,nNCases] ...
    = size(r1Results);

averagedR1 = squeeze(mean(r1Results,3));

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
nearestNeighbours = getValuesFromStringEnumeration( ...
    configuration.nearestNeighbourCases,';','numeric')'

reproductionResults(1,:) = nearestNeighbours;
reproductionResults(2,:) = averagedRates;
reproductionResults(3,:) = squeeze(mean(mean(averageDifferences,2),1));
reproductionResults(4,:) = averageEstimatedR1;

material = strsplit(fileName,'_');
material = material{1};

figure(1)
hold on
legendEntries = {};
for referenceOrientationNr = 1:orientationsCount
    for referencePositionNr = 1:positionsCount
        plot(nearestNeighbours,squeeze(averageDifferences( ...
            referenceOrientationNr,referencePositionNr,:)) ...
            ,'*-','LineWidth',1.5)
        orientationAngle = rad2deg(orientationAngles( ...
            referenceOrientationNr));
        positionAngle = rad2deg(positionAngles( ...
            referencePositionNr));
        legendEntries{end+1} = sprintf("$\\Theta: %.1f$ $\\varphi: %.1f$" ...
            ,orientationAngle,positionAngle); %#ok<SAGROW>
    end
end
grid minor
legend(legendEntries,'Interpreter','Latex')
multiLineTitle = {sprintf('%s Reproduce R_1 with less nearest',material) ...
    , sprintf(' neighbours (Data: %s)' ,creationDate)};
title(multiLineTitle)
xlabel('Number of nearest neighbours to repoduce R_1')
ylabel('Averaged difference to original data')
creationDate = datestr(now,'yyyymmdd');
filePath = sprintf('%s%s_%s_%s.png',path2Results,creationDate ...
    ,material,'ReproduceRelaxationRates');
saveas(gcf,filePath);


