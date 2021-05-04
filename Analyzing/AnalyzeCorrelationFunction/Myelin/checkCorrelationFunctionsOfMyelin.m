
%% load data
clc
clear all
close all

path2Load = ['C:\Users\maxoe\Google Drive\Promotion\Results\Myelin' ...
    '\Relaxation\Lipid_H_500ns_1ps_wh_6932_resultsOriDep_20210430.mat'];

data = load(path2Load);
B0 = data.B0;

path2Save = ['C:\Users\maxoe\Google Drive\Promotion\Results\Myelin' ...
    '\Correlation\CorrelationFunctionsAt' num2str(B0) 'T.fig'];

sumCorrelationFunctionW0 = data.correlationFunctionW0Saver;
sumCorrelationFunction2W0 = data.correlationFunction2W0Saver;

relaxationRates = data.r1WithPerturbationTheory;
[orientationCount,positionCount,~] = size(relaxationRates);
try
    atomsWithCalculatedR1Rates = (relaxationRates > 0.000001);
    atomsWithCalculatedR1Rates = logical(squeeze(mean(mean( ...
        atomsWithCalculatedR1Rates,1),2))');
    atomCount = data.atomCounter;
catch
    atomCount = sum(atomsWithCalculatedR1Rates);
end


deltaT = data.deltaT;

%% analyze data
correlationFunctionW0 = sumCorrelationFunctionW0/atomCount;
correlationFunction2W0 = sumCorrelationFunction2W0/atomCount;

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
figs(1) = figure(1)
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
figs(2) = figure(2)
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

savefig(figs,path2Save)




