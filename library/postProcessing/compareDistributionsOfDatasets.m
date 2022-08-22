function compareDistributionsOfDatasets(dataset_path,saving)

dataset = load(dataset_path);
binWidth = 0.1;
startDateOfSimulation = dataset.startDateOfSimulation;
whichLipid = dataset.whichLipid;
oriAngles = rad2deg(dataset.orientationAngles);
posAngles = rad2deg(dataset.positionAngles);
atomIndex = dataset.atomIndex;
numberOfAtoms = size(dataset.meanPositions,1);
calculatedAtomLocations = dataset.meanPositions(atomIndex,:);
initializeFigure('legend',false);
histogram(calculatedAtomLocations(:,1),6);
title(sprintf('%s: validation dataset (%i atoms), $\\theta$: %.2f $\\varphi$: %.2f', ...
    whichLipid,size(calculatedAtomLocations,1),oriAngles(1),posAngles(1)));
ylabel('Frequency');
xlabel('x-location [nm]');
if saving
    savingPath = initializeSystemForSavingDistributionPlots();
    savingName = sprintf('%s_DistributuionOfAtomsAlongX_validationDataset_%s' ...
        ,startDateOfSimulation,whichLipid);
    print(gcf,[savingPath savingName],'-dpng','-r300');
end

initializeFigure('legend',false);
title(sprintf('%s: whole dataset (%i atoms), $\\theta$: %.2f $\\varphi$: %.2f', ...
    whichLipid,numberOfAtoms,oriAngles(1),posAngles(1)));
histogram(dataset.meanPositions(:,1),6);
ylabel('Frequency');
xlabel('x-location [nm]');
if saving
    savingPath = initializeSystemForSavingDistributionPlots();
    savingName = sprintf('%s_DistributuionOfAtomsAlongX_wholeDataset_%s' ...
        ,startDateOfSimulation,whichLipid);
    print(gcf,[savingPath savingName],'-dpng','-r300');
end




end
