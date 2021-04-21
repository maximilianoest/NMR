
%% load data
dataLoaded = true;
if ~exist('correlationFunctionW0','var') || dataLoaded
    clc
    clear all
    close all
    loadedData = load(['C:\Users\maxoe\Documents\ServerData' ...
        '\Lipid_ExampleTrajSet_nH676_resultsOrientationDependency.mat']);
    correlationFunctionW0 = loadedData.correlationFunctionW0Saver;
    correlationFunction2W0 = loadedData.correlationFunction2W0Saver;
end

%% analyze data

averageCorrelationFunctionW0 = squeeze(mean(correlationFunctionW0,3));
averageCorrelationFunction2W0 = squeeze(mean(correlationFunction2W0,3));

effectiveCorrelationFunctionW0 = squeeze(mean( ...
    averageCorrelationFunctionW0,2));
effectiveCorrelationFunction2W0 = squeeze(mean( ...
    averageCorrelationFunction2W0,2));

%% plotting 
orientationsCount = size(effectiveCorrelationFunctionW0,1);
orientationAngles = linspace(0,90,orientationsCount);
tauAxis = 0:1e-12:(size(averageCorrelationFunctionW0,3)-1)*1e-12;
legendEntries = {};
figure(1)
hold on
for i = 1:orientationsCount
    plot(tauAxis,abs(real(effectiveCorrelationFunctionW0(i,:) ...
        /effectiveCorrelationFunctionW0(i,1))))
    legendEntries{i} = num2str( ...
        orientationAngles(i),'Theta: %.2f'); %#ok<SAGROW>
    axis([0 1.4e-7 0 0.5])
    pause(0.2)
end
hold off
grid on
legend(legendEntries)
title('Efffective Correlation Function at w0')
xlabel('tau')

orientationsCount = size(effectiveCorrelationFunction2W0,1);
legendEntries = {};
figure(2)
hold on
for i = 1:orientationsCount
    plot(tauAxis,abs(real(effectiveCorrelationFunction2W0(i,:) ...
        /effectiveCorrelationFunction2W0(i,1))))
    legendEntries{i} = num2str( ...
        orientationAngles(i),'Theta: %.2f'); %#ok<SAGROW>
    axis([0 1.4e-7 0 0.5])
    pause(0.2)
end
hold off
grid on
legend(legendEntries)
title('Effective Correlation Function at w0')
xlabel('tau')

%% save all variables except for correlationFunctions

fieldNames = fieldnames(loadedData);
% dataWithoutCorrelationFunctions = struct([]);

for fieldNr = 1:length(fieldNames)
    
    switch fieldNames{fieldNr}
          case {'correlationFunctionW0Saver' ...
                  ,'correlationFunction2W0Saver'}
          otherwise
            dataWithoutCorrelationFunctions(1).(fieldNames{fieldNr}) = ...
                loadedData.(fieldNames{fieldNr});
    end
end

save(['C:\Users\maxoe\Google Drive\Promotion\Data\Myelin' ...'
    '\Lipid_RelaxationRates_nH40.mat'],'dataWithoutCorrelationFunctions');



