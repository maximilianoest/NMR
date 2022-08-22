clear all; clc; close all;

%% first
% close all
colorArray = [252, 3, 3;252, 194, 3;82, 252, 3;3, 252, 240;3, 252, 240; ...
    173, 3, 252;252, 3, 181;3, 186, 252;231, 252, 3;3, 152, 252]/252;

baseConfiguration = readConfigurationFile('../baseConfiguration.txt');
addpath(genpath(baseConfiguration.path2LibraryOnLocalMachine));
results = load("C:\Users\maxoe\Google Drive\Promotion\Results\LIPIDS\PSM\orientationDependency\20220530_Results_orientationDependency_PSMlipid.mat");

saving = 1;
if saving
    savingPath = initializeSystemForSavingPlots('locationDependency' ...
    ,results.whichLipid);
end
atomIndex = results.atomIndex(1:results.atomCounter);
positionAnglesIndex = 3;
orientationAngleIndex = 1;
positionAngle = results.positionAngles(positionAnglesIndex)
orientationAngle = results.orientationAngles(orientationAngleIndex)
zAxis = [0 0 1];
rotationMatrixPosition = get3DRotationMatrix( ...
    positionAngle,zAxis);
yAxis = [0 1 0];
rotationMatrixOrientation = get3DRotationMatrix( ...
    orientationAngle,yAxis);
totalRotationMatrix = ...
    rotationMatrixOrientation*rotationMatrixPosition;
untouchedMeanPositions = results.meanPositions(atomIndex,:);
meanPositions = (totalRotationMatrix*untouchedMeanPositions')';

allR1 = results.scaledR1Rates;
allR1 = allR1(:,:,1:results.atomCounter);
% overallR1 = squeeze(mean(mean(allR1,2),1))';
overallR1 = squeeze(allR1(orientationAngleIndex,positionAnglesIndex,:));

dimensions = ["X", "Y", "Z"];

initializeFigure('legend',false);
plot3(meanPositions(:,1),meanPositions(:,2),meanPositions(:,3),'*');
xlabel('X');
ylabel('Y');
zlabel('Z');

for dimension = 1:size(dimensions,2)
    initializeFigure('legend',false);
    scatter(meanPositions(:,dimension),overallR1,'LineWidth',1.3);
    title(sprintf('Relaxation rate R$_1$ in dimension %s (%s)' ...
        ,dimensions(dimension),results.whichLipid));
    xlabel(sprintf('Location in %s direction',dimensions(dimension)));
    if saving
        savingName = sprintf('scatterPlotLocationDepR1_Dim%s_%s' ...
            ,dimensions(dimension),results.whichLipid);
        print(gcf,[savingPath savingName],'-dpng','-r300');
    end
end

% THE RELAXATION RATES AT THE TAILS ARE SMALLER THAN AT THE HEADS LIKE
% PRESENTED BY SCHYBOLL 2019
% Comparison with Fig. 20.16 in Levitt: Heads will be positioned more to
% the minimum of the curve while the tails will be placed more to lower
% correlation times what is in accordance with the results of this
% investigation.

%% second
% averaged over regions:
initializeFigure();
legendEntries = {};
locationSteps = 25;

averageR1 = zeros(1,locationSteps);
averagePositions = zeros(1,locationSteps);

for dimension = 1:size(dimensions,2)
    positionsHydrogenDim = meanPositions(:,dimension);
    minLocation = floor(min(positionsHydrogenDim));
    maxLocation = ceil(max(positionsHydrogenDim));
    locationsToAverage = linspace(minLocation,maxLocation,locationSteps);
    for locationNr = 1:size(locationsToAverage,2)-1
        indices = (positionsHydrogenDim>=locationsToAverage(locationNr)) & ...
            (positionsHydrogenDim<=locationsToAverage(locationNr+1));
        averageR1(locationNr) = mean(overallR1(indices));
    end
    plot(locationsToAverage,averageR1);
    legendEntries{end+1} = sprintf('Dimension: %s',dimensions(dimension)); %#ok<SAGROW>
end
legend(legendEntries);

title(sprintf('Relaxation rate dependent on location in lipid (%s)' ...
    ,results.whichLipid));
if saving
    savingName = sprintf('locationDependentR1_%s',results.whichLipid);
    print(gcf,[savingPath savingName],'-dpng','-r300');
end


%% YOU CAN DETERMINE R1 OF HEADS IF KNOWN POSITIONS
% leftEnd = 2.2;
% rightEnd = 5.9;
% headPositions = (meanPositions(:,1) <= 2.2 | meanPositions(:,1) >= 5.9);
% r1OfHead = overallR1(headPositions);
% for dimension = 1:1 %size(dimensions,2)
%     initializeFigure('legend',false);
%     scatter(meanPositions(headPositions,dimension) ...
%         ,r1OfHead,'LineWidth',1.3);
%     title(sprintf('Relaxation rate R$_1$ of heads in dimension %s (%s)' ...
%         ,dimensions(dimension),results.whichLipid));
%     xlabel(sprintf('Location in %s direction',dimensions(dimension)));
%     if saving
%         savingName = sprintf('scatterPlotLocationDepR1ofHeads_Dim%s_%s' ...
%             ,dimensions(dimension),results.whichLipid);
%         print(gcf,[savingPath savingName],'-dpng','-r300');
%     end
% end
% 
% initializeFigure('legend',false);
% histogram(r1OfHead,'BinWidth',1);
% xlabel('Relaxation rate $R_1$ [Hz]');
% ylabel('Frequency');
% title('Distribution of relaxation rates at the lipid heads');
% if saving
%     savingName = sprintf('histogramR1LipidHeads_Dim%s_%s' ...
%         ,'X',results.whichLipid);
%     print(gcf,[savingPath savingName],'-dpng','-r300');
% end


% You see there is nearly no change of R1 in Y and Z direction 





% diffusion coefficient

% minus average to get the same zentral position for each lipid