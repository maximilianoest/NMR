
%% load data
clc
clear all
close all
loadedData = load(['C:\Users\maxoe\Google Drive\Promotion' ...
    '\Results\Myelin' ...
    '\Lipid_H_500ns_1ps_nH1000_180Hs_resultsOriDep_20210426.mat']);
correlationFunctionW0 = loadedData.correlationFunctionW0Saver;
correlationFunction2W0 = loadedData.correlationFunction2W0Saver;

%% analyze data

relaxationRates = loadedData.averageRelaxationRates;
effectiveRelaxationRates = squeeze(mean(relaxationRates,2));

precision = realmin('double');
numberOfAtoms = sum(effectiveRelaxationRates(1,:) > 100*precision);

effectiveCorrelationFunctionW0 = squeeze(mean( ...
    correlationFunctionW0,2));
effectiveCorrelationFunction2W0 = squeeze(mean( ...
    correlationFunction2W0,2));

save('correlationFunctionsMyelin_nH164' ...
    ,'effectiveCorrelationFunctionW0' ...
    ,'effectiveCorrelationFunction2W0','-v7.3');

%% plotting 
orientationsCount = size(effectiveCorrelationFunctionW0,1);
orientationAngles = linspace(0,90,orientationsCount);
tauAxis = 0:2*1e-12:(size(correlationFunctionW0,3)-1)*2*1e-12;
legendEntries = {};
fig(1) = figure(1)
hold on
for i = 1:orientationsCount
    plot(tauAxis,abs(real(effectiveCorrelationFunctionW0(i,:) ...
        /effectiveCorrelationFunctionW0(i,1))))
    legendEntries{i} = num2str( ...
        orientationAngles(i),'Theta: %.2f'); %#ok<SAGROW>
    axis([0 2.8e-7 0 0.5])
    pause(0.2)
end
hold off
grid on
legend(legendEntries)
title('Effective Correlation Function at w0')
xlabel('tau')

orientationsCount = size(effectiveCorrelationFunction2W0,1);
legendEntries = {};
fig(2) = figure(2)
hold on
for i = 1:orientationsCount
    plot(tauAxis,abs(real(effectiveCorrelationFunction2W0(i,:) ...
        /effectiveCorrelationFunction2W0(i,1))))
    legendEntries{i} = num2str( ...
        orientationAngles(i),'Theta: %.2f'); %#ok<SAGROW>
    axis([0 2.8e-7 0 0.5])
    pause(0.2)
end
hold off
grid on
legend(legendEntries)
title('Effective Correlation Function at 2w0')
xlabel('tau')

savefig(fig,'CorrelationFunctions.fig')




