
%% load data
dataLoaded = true;
if ~exist('correlationFunctionW0','var') || dataLoaded
    clc
    clear all
    close all
    savedData = load(['C:\Users\maxoe\Documents\ServerData' ...
        '\Lipid_ExampleTrajSet_nH676_resultsOrientationDependency.mat']);
    correlationFunctionW0 = savedData.correlationFunctionW0Saver;
    correlationFunction2W0 = savedData.correlationFunction2W0Saver;
    clearvars -except correlationFunctionW0 correlationFunction2W0
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



