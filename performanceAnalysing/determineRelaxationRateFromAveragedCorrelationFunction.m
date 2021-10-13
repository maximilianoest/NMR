% clc
% clear all


results = load("C:\Users\maxoe\Google Drive\Promotion\Results\nearestNeighboursAnalysis\Server\20211001_Results_nearestNeighbourAnalysis_Lipid_H_500ns_4ps_wh.mat");
results = load("C:\Users\maxoe\Google Drive\Promotion\Results\Myelin\Relaxation\035Tesla\20210522_Lipid_H_500ns_4ps_wh_resultsOriDep.mat");

effectiveCorrelationFunction1W0 = squeeze(mean(mean(results.correlationFunction1W0Saver,3),2));
effectiveCorrelationFunction1W0 = effectiveCorrelationFunction1W0./effectiveCorrelationFunction1W0(:,1);
effectiveCorrelationFunction2W0 = squeeze(mean(mean(results.correlationFunction2W0Saver,3),2));
effectiveCorrelationFunction2W0 = effectiveCorrelationFunction2W0./effectiveCorrelationFunction2W0(:,1);

figure(1)
plot(abs(effectiveCorrelationFunction1W0(1,:)),'LineWidth',1.5)
hold on
plot(abs(effectiveCorrelationFunction2W0(1,:)),'LineWidth',1.5)
plot(abs(effectiveCorrelationFunction1W0(2,:)),'LineWidth',1.5)
plot(abs(effectiveCorrelationFunction2W0(2,:)),'LineWidth',1.5)
hold off
legend('\omega_0 0°','2\omega_0 0°','\omega_0 90°','2\omega_0 90°')
grid minor
title('Korrelationsfunktion Mylein')
