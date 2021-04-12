clc 
clear all
close all


data = load('C:\Users\maxoe\Google Drive\Promotion\Results\PLPC\PLPC_allData_CorrelationFunction.mat');
meanData = data.meanCorrelation;
normData = data.normCorrelation;


figure('Name','Normalized Correlation Functions PLPC')
timeAxis = 0:1e-12:(size(meanData,2)-1)*1e-12;
plot(timeAxis,normData(1,:),'LineWidth',1.2)
hold on
plot(timeAxis,normData(2,:),'LineWidth',1.2)
hold off
xlabel('Tau [s]')
grid on
axis([min(timeAxis) max(timeAxis) 0 0.5])
title('Normalized Correlation Functions PLCP')