
%% load data
dataLoaded = true;
if ~exist('correlationFunctionW0','var')
    clc
    clear all
    close all
    savedData = load(['C:\Users\maxoe\Google Drive\Promotion\Results\' ...
        'Myelin\'...
        'Lipid_ExampleTrajSet_nH676_resultsOrientationDependency.mat']);
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
figure(1)
hold on
for i = 1:size(effectiveCorrelationFunctionW0,1)
    plot(abs(real(effectiveCorrelationFunctionW0(i,:))))
    pause(0.2)
end
hold off
grid on
title('Efffective Correlation Function at w0')
xlabel('tau')




