clc; close all;

% Questions:
% - how are the indices organized in a dataset
% 

baseConfiguration = readConfigurationFile('../baseConfiguration.txt');
addpath(genpath(baseConfiguration.path2LibraryOnLocalMachine));

dataLoaded = exist('fileList','var');
if ~dataLoaded
    dataFolder = 'C:\Users\maxoe\Google Drive\Promotion\Data\Lipids\DOPS\singleLipids\';
    fileList = dir(sprintf('%s*.mat',dataFolder));
    for fileNr = 1:length(fileList)
        fileList(fileNr).data = load(sprintf('%s%s',dataFolder ...
            ,fileList(fileNr).name));
    end
end

saving = 1;
savingPath = initializeSystemForSavingPlots('hydrogenMovement' ...
    ,results.whichLipid);

dt = 1;
timeSteps = 1:dt:1;%round(size(fileList(1).data.trajectories,3)*0.2);
fig = initializeFigure('legend',false,'lineWidth',1);
view(-37.5,30);
axis([0 7 0 7 0 7]);
colorArray = [252, 3, 3;252, 194, 3;82, 252, 3;3, 252, 240;3, 252, 240; ...
    173, 3, 252;252, 3, 181;3, 186, 252;231, 252, 3;3, 152, 252]/252;

plt = {};
% the first 8 atoms can be identified as hydrogen atoms from the head due
% to structural formula in masters thesis
sepIdx=8;
for timeStep = timeSteps
    for lipidNr = 1:length(fileList)
        
        x = squeeze(fileList(lipidNr).data.trajectories(:,1,timeStep));
        y = squeeze(fileList(lipidNr).data.trajectories(:,2,timeStep));
        z = squeeze(fileList(lipidNr).data.trajectories(:,3,timeStep));
        plt.(sprintf('headPlot%i',lipidNr)) = plot3(x(1:sepIdx),y(1:sepIdx) ...
            ,z(1:sepIdx),'o','Color',colorArray(lipidNr,:));
        plt.(sprintf('tailPlot%i',lipidNr)) = plot3(x(sepIdx+1:end) ...
            ,y(sepIdx+1:end),z(sepIdx+1:end),'*' ...
            ,'Color',colorArray(lipidNr,:));
        [caz,cel] = view;
        view([caz,cel]);
        xlabel('X');
        ylabel('Y');
        zlabel('Z');
    end
    shg;
    pause(0.2);
    plotNames = fieldnames(plt);
    for lipidNr = 1:length(plotNames)
        if timeStep ~= timeSteps(end)
            delete(plt.(plotNames{lipidNr}));
        end
    end
    
end

% CONSEQUENCE: There is defenitly more movements within the tails than
% within the heads of the lipids.
% interpretation: relaxation rate at the head should smaller than at the
% tails, what is different to Shyboll et al. 2019
% take a look at the dataset from last calculations and look at the
% location dependent R1

results = load("C:\Users\maxoe\Google Drive\Promotion\Results\LIPIDS\DOPS\orientationDependency\20220307_Results_orientationDependency_DOPSlipid.mat");
r1 = results.scaledR1Rates;
r1 = r1(:,:,1:results.atomCounter);
overallR1 = squeeze(mean(mean(r1,2),1))';
meanPositions = results.meanPositions;
atomIndex = results.atomIndex(1:results.atomCounter);
meanPositions = meanPositions(atomIndex,:);

dimensions = ["X", "Y", "Z"];

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

% You see there is nearly no change of R1 in Y and Z direction 





% diffusion coefficient

% minus average to get the same zentral position for each lipid