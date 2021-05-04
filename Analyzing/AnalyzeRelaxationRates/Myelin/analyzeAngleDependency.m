clc
clear all %#ok<CLALL>

path2Load = ['C:\Users\maxoe\Google Drive\Promotion\Results\Myelin' ...
    '\Relaxation\Lipid_H_500ns_1ps_wh_6932_resultsOriDep_20210430.mat'];

data = load(path2Load);
B0 = data.B0;

path2Save = ['C:\Users\maxoe\Google Drive\Promotion\Results\Myelin' ...
    '\Relaxation\RelaxationRatesAt' num2str(B0) 'T.fig'];

r1Perturbation = data.r1WithPerturbationTheory;
[orientationCount,positionCount,~] = size(r1Perturbation);
try
    atomsWithCalculatedR1Rates = (r1Perturbation > 0.000001);
    atomsWithCalculatedR1Rates = logical(squeeze(mean(mean( ...
        atomsWithCalculatedR1Rates,1),2))');
    atomCount = data.atomCounter;
catch
    atomCount = sum(atomsWithCalculatedR1Rates);
end

calculatedR1Rates = r1Perturbation(:,:,1:atomCount);
averageRelaxationRates = mean(calculatedR1Rates,3);
medianRelaxationRates = median(calculatedR1Rates,3);

effectiveRelaxationRatesMean = mean(averageRelaxationRates,2);
tmp = reshape(calculatedR1Rates,orientationCount,atomCount*positionCount);
effectiveRelaxationRatesMedian = median(tmp,2);

figs(1) = figure(1);
orientationAngles = linspace(0,90,orientationCount);
plot(orientationAngles,effectiveRelaxationRatesMean,'LineWidth',1.5);
hold on
plot(orientationAngles,effectiveRelaxationRatesMedian,'LineWidth',1.5);
hold off
legend('Mean','Median')
grid on
title(['Angle Dependency of Relaxation Rate at ' num2str(B0) ' Tesla'])
xlabel('Angle [°]')
ylabel('Relaxation Rates [Hz]')

savefig(figs,path2Save)


