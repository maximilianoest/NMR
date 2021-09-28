clc
clear all  %#ok<CLALL>
close all

results = load("C:\Users\maxoe\Google Drive\Promotion\Results\OrientationDependency\20210928_Results_OrientationDependency_Lipid_H_500ns_1ps_nH40.mat");

relaxationRates = results.r1WithPerturbationTheory; 
atomCounter = results.atomCounter;

averageR1 = squeeze(mean(relaxationRates(:,:,1:atomCounter,:),3));
effectiveR1 = squeeze(mean(averageR1,2));
overallR1 = squeeze(mean(effectiveR1,1));

nearestNeighbourCases = results.nearestNeighbourCases;

figure(1)
plot(nearestNeighbourCases,overallR1,'LineWidth', 1.5)
xlabel('Nearest Neighbours')
ylabel('Overall R1 [Hz]')
title('Necessity of number of nearest Neighbours')
grid minor



