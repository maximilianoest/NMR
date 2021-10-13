clc
clear all
close all

%% Analyse whether the same results can be achieved with single precision
configuration = readConfigurationFile('..\config.txt');
path2Conclusions = [configuration.path2ResultsOnLocalMachine ...
     configuration.kindOfResults '\Conclusions\'];

path2OldData = "C:\Users\maxoe\Google Drive\Promotion\Results\nearestNeighboursAnalysis\Server\20211001_Results_nearestNeighbourAnalysis_Lipid_H_500ns_4ps_wh.mat";
oldData = load(path2OldData);
oldR1 = oldData.r1WithPerturbationTheory;
oldAtomCounter = oldData.atomCounter;

path2NewData = "C:\Users\maxoe\Google Drive\Promotion\Results\PerformanceAnalysis\Server\20211008_Results_performanceAnalysing_Lipid_H_500ns_4ps_wh.mat";
newData = load(path2NewData);
newR1 = newData.r1WithPerturbationTheory;
newTimeTracks = newData.timeTracks;
newAtomCounter = newData.atomCounter;

whichNearestNeighbourCase = find( ...
    oldData.nearestNeighbourCases == newData.nearestNeighbours);
oldR1 = oldR1(:,:,1:newAtomCounter,whichNearestNeighbourCase);
oldAveragedR1 = mean(oldR1,3);
oldEffectiveR1 = squeeze(mean(oldAveragedR1,2));

newAveragedR1 = mean(newR1,3);
newEffectiveR1 = mean(newAveragedR1,2);
orientationAngles = rad2deg(newData.orientationAngles);
dateNow = datestr(now,'yyyymmdd_HHMM_');

figure(1)
hold on
plot(orientationAngles,oldAveragedR1(:,1),'--*')
plot(orientationAngles,oldAveragedR1(:,2),'--*')
plot(orientationAngles,newAveragedR1(:,1),'-.o')
plot(orientationAngles,newAveragedR1(:,2),'-.o')
hold off
legend('double precision \phi = 0°','double precision \phi = 90°' ...
    ,'single precision \phi = 0°','single precision \phi = 90°' ...
    ,'Location','NorthWest')
legend boxoff
title(['Averaged relaxation rates with single and double precision' ...
    ' for 1000 atoms.'])
grid minor
saveas(gcf,[path2Conclusions dateNow ...
    'ComparisonSingleDoublePrecision.png'])
save([path2Conclusions dateNow 'Data_Comparison.mat'])
close all





