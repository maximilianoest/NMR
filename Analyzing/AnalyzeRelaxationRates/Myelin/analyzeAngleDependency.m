clc
clear all %#ok<CLALL>

data = load(['C:\Users\maxoe\Google Drive\Promotion\Results\Myelin' ...
    '\Relaxation\Water_H_50ns_025ps_nH2500_resultsOriDep_20210427.mat']);

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

calculatedR1Rates = r1Perturbation(:,:,atomsWithCalculatedR1Rates);
averageRelaxationRates = mean(calculatedR1Rates,3);
medianRelaxationRates = median(calculatedR1Rates,3);

effectiveRelaxationRatesMean = mean(averageRelaxationRates,2);
tmp = reshape(calculatedR1Rates,orientationCount,atomCount*positionCount);
effectiveRelaxationRatesMedian = median(tmp,2);

figure(1)
orientationAngles = linspace(0,90,orientationCount);
plot(orientationAngles,effectiveRelaxationRatesMean);
hold on
plot(orientationAngles,effectiveRelaxationRatesMedian);
hold off
legend('Mean','Median')


