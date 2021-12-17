clc
clear all
%#ok<*SAGROW>

results = load("C:\Users\maxoe\Google Drive\Promotion\Results\determineSamplingRate\Server\20211027_Results_determineSamplingRate_Lipid.mat");

r1WithExplIntegration = results.r1WithExplIntegration;
dataSets = fieldnames(r1WithExplIntegration);
calculatedAtoms = results.atomsToCalculate;

for dataSetNumber = 1:length(dataSets)
    dataSetNames(dataSetNumber,:) = string(dataSets{dataSetNumber});
    r1Results(dataSetNumber,:,:,:) = ...
        r1WithExplIntegration.(dataSetNames(dataSetNumber,:)).results;
    averageR1Results(dataSetNumber,:,:) = ...
        mean(r1Results(dataSetNumber,:,:,:),4);
    effectiveR1Results(dataSetNumber,:) = ...
        mean(averageR1Results(dataSetNumber,:,:),3);
end

orientationAngles = rad2deg(results.orientationAngles);

figure(1)
hold on
legendEntries = {};
for dataSetNumber = 1:length(dataSets)
    dataSetName = dataSetNames(dataSetNumber);
    legendEntries{dataSetNumber} = dataSetName;
    plot(orientationAngles,effectiveR1Results( ...
        dataSetNumber,:),'LineWidth',1.5);
end
hold off
grid minor
legend(legendEntries,'Interpreter','none','Location','northwest')
title('Myelin orientation dependency for different sampling rates')
xlabel('Orientation angle [°]')
ylabel('Relaxation rate [Hz]')


