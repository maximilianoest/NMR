
%% load data
clc
clear all
close all
loadedData = load(['C:\Users\maxoe\Google Drive\Promotion\Results' ...
    '\Myelin\Relaxation' ...
    '\Lipid_H_500ns_1ps_nH1000_resultsOriDep_20210426.mat']);
sumCorrelationFunctionW0 = loadedData.correlationFunctionW0Saver;
sumCorrelationFunction2W0 = loadedData.correlationFunction2W0Saver;
relaxationRates = loadedData.averageRelaxationRates;
try
    atomsCount = loadedData.atomCounter;
catch
    warning('Atom Count not defined')
    calculatedRelaxationRates = squeeze(relaxationRates(1,1,:));
    relaxationRatesNotZero = calculatedRelaxationRates > 0.0001;
    atomsCount = sum(relaxationRatesNotZero);
end


deltaT = loadedData.deltaT;

%% analyze data
correlationFunctionW0 = sumCorrelationFunctionW0/atomsCount;
correlationFunction2W0 = sumCorrelationFunction2W0/atomsCount;

effectiveRelaxationRates = squeeze(mean(relaxationRates,2));

effectiveCorrelationFunctionW0 = squeeze(mean( ...
    correlationFunctionW0,2));
effectiveCorrelationFunction2W0 = squeeze(mean( ...
    correlationFunction2W0,2));

%% plotting 
orientationsCount = size(effectiveCorrelationFunctionW0,1);
orientationAngles = linspace(0,90,orientationsCount);
tauAxis = 0:deltaT:(size(sumCorrelationFunctionW0,3)-1)*deltaT;
tauMin = 0;
tauMax = 5e-7;
valueMin = 0;
valueMax = 0.5;
legendEntries = {};
fig(1) = figure(1)
hold on
for i = 1:orientationsCount
    plot(tauAxis,abs(real(effectiveCorrelationFunctionW0(i,:) ...
        /effectiveCorrelationFunctionW0(i,1))),'LineWidth',1.5)
    legendEntries{i} = num2str( ...
        orientationAngles(i),'Theta: %.2f'); %#ok<SAGROW>
    axis([tauMin tauMax valueMin valueMax])
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
        /effectiveCorrelationFunction2W0(i,1))),'LineWidth',1.5)
    legendEntries{i} = num2str( ...
        orientationAngles(i),'Theta: %.2f'); %#ok<SAGROW>
    axis([tauMin tauMax valueMin valueMax])
end
hold off
grid on
legend(legendEntries)
title('Effective Correlation Function at 2w0')
xlabel('tau')

savefig(fig,'CorrelationFunctions.fig')




