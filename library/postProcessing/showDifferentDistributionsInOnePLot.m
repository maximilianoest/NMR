function showDifferentDistributionsInOnePLot(validationSet_path, ...
    scalFacSet_path,dimension,binWidth,saving)

switch dimension
    case 'X'
        indxDimension = 1;
    case 'Y'
        indxDimension = 2;
    case 'Z'
        indxDimension = 3;
    otherwise
        fprintf('Dimension %s not found!!',dimension);
        error('dimensionNotFound:showDifferentDistributionsInOnePlot', ...
            'The dimension %s cannot be found',dimension);
end

% compare the datasets
validationDataset = load(validationSet_path);

startDateOfSimulation = validationDataset.startDateOfSimulation;
whichLipid = validationDataset.whichLipid;
fprintf('Lipid of validation set: %s\n',whichLipid);
oriAngles = rad2deg(validationDataset.orientationAngles);
posAngles = rad2deg(validationDataset.positionAngles);
validationAtomIndex = validationDataset.atomIndex;

validationCalculatedAtomLocations = validationDataset.meanPositions(validationAtomIndex,indxDimension);
wholeSetAtomLocations = validationDataset.meanPositions(:,indxDimension);

scalFacDataset = load(scalFacSet_path);
whichLipid = scalFacDataset.whichLipid;
fprintf('Lipid of scaling factors set: %s\n',whichLipid);
scalFacCalculatedAtomLocations = scalFacDataset.meanPositions( ...
    scalFacDataset.atomIndex,indxDimension);
initializeFigure();
legendEntries = {};
hold on
h1 = histogram(validationCalculatedAtomLocations);
h1.Normalization = 'probability';
h1.LineStyle = '-';
h1.BinWidth = binWidth;
legendEntries{end+1} = sprintf('%i atoms',size(validationCalculatedAtomLocations,1));
h2 = histogram(wholeSetAtomLocations);
h2.Normalization = 'probability';
h2.LineStyle = '--';
h2.BinWidth = binWidth;
legendEntries{end+1} = sprintf('%i atoms',size(wholeSetAtomLocations,1));
hold off
legend(legendEntries);
title(sprintf('%s distribution of H atoms, $\\theta$: %.2f $\\varphi$: %.2f', ...
    whichLipid,oriAngles(1),posAngles(1)));
ylabel('Relative frequency');
xlabel('x-location [nm] head $\rightarrow$ tail tail $\leftarrow$ head');

if saving
    savingPath = initializeSystemForSavingDistributionPlots();
    savingName = sprintf('%s_HistoDistributionOfAtomsComparison_%s' ...
        ,startDateOfSimulation,whichLipid);
    print(gcf,[savingPath savingName],'-dpng','-r300');
end

figure('visible','off');
h3 = histogram(scalFacCalculatedAtomLocations);
h3.Normalization = 'probability';
h3.BinWidth = binWidth;



validationHistoData = h1.Values;
ValidationSetPositionAxis = h1.BinEdges(1:end-1);
wholeSetHistoData = h2.Values;
wholeSetPositionAxis = h2.BinEdges(1:end-1);
scalFacHistoData = h3.Values;
scalFacPositionAxis = h3.BinEdges(1:end-1);


initializeFigure();
legendEntries = {};
hold on
plot(ValidationSetPositionAxis,validationHistoData,'o--');
legendEntries{end+1} = sprintf('%i atoms',size( ...
    validationCalculatedAtomLocations,1));
plot(wholeSetPositionAxis,wholeSetHistoData,'o--');
legendEntries{end+1} = sprintf('%i atoms',size( ... 
    wholeSetAtomLocations,1));
plot(scalFacPositionAxis,scalFacHistoData,'o--');
legendEntries{end+1} = sprintf('%i atoms',size( ...
    scalFacCalculatedAtomLocations,1));
hold off
legend(legendEntries)
title(sprintf('%s distribution of H atoms, $\\theta$: %.2f $\\varphi$: %.2f', ...
    whichLipid,oriAngles(1),posAngles(1)));
ylabel('Relative frequency');
xlabel('x-location [nm] head $\rightarrow$ tail tail $\leftarrow$ head');
if saving
    savingPath = initializeSystemForSavingDistributionPlots();
    savingName = sprintf('%s_LineDistributionOfAtomsComparison_%s' ...
        ,startDateOfSimulation,whichLipid);
    print(gcf,[savingPath savingName],'-dpng','-r300');
end


% smallestLocation = min(wholeSetAtomLocations(:,indxDimension));
% largestLocation = max(wholeSetAtomLocations(:,indxDimension));
% binEdges = smallestLocation;
% largestLocationReached = false;
% % determine edges of bins
% while ~largestLocationReached
%     binEdges(end+1) = binEdges(end) + binWidth;
%     if binEdges(end) > largestLocation 
%         largestLocationReached = true;
%         binEdges(end) = binEdges(end)+0.0001;
%     end
% end
% 
% % find data for bins
% for binEdgeCount = 1:length(binEdges)-1
%     lowerEdge = binEdges(binEdgeCount);
%     upperEdge = binEdges(binEdgeCount+1);
%     relativeNumberOfAtomsInBinWholeSet(binEdgeCount) = ...
%         getRelativeNumberOfAtomsInBin(lowerEdge,upperEdge, ...
%         wholeSetAtomLocations);
%     relativeNumberOfAtomsInBinValidationSet(binEdgeCount) = ...
%         getRelativeNumberOfAtomsInBin(lowerEdge,upperEdge, ...
%         validationCalculatedAtomLocations);
%     relativeNumberOfAtomsInBinScalFacSet(binEdgeCount) = ...
%         getRelativeNumberOfAtomsInBin(lowerEdge,upperEdge, ...
%         scalFacCalculatedAtomLocations);
% end



end
