clc
clear all
close all

data = load('C:\Users\maxoe\Google Drive\Promotion\Results\PSM\PSM_LIPIDS_H_250ns_1ps_resultsRelaxationRates.mat');
atomNumber = 160;

r1WithLipariSzabo = data.r1WithLipariSzabo;
r1WithPerturbationTheory = data.r1WithPerturbationTheory;
r1WithSchroedingerEquation = data.r1WithSchroedingerEquation;

meanR1WithPerturbationTheory = data.meanR1WithPerturbationTheory;
meanR1WithLipariSzabo = data.meanR1WithLipariSzabo;
meanR1WithSchroedingerEquation = data.meanR1WithSchroedingerEquation;

figure(1)
plot(r1WithPerturbationTheory(1:atomNumber),'b','LineWidth',1.5)
hold on
plot(meanR1WithPerturbationTheory(1:atomNumber),'--b' ...
    ,'LineWidth',1.5)
plot(r1WithLipariSzabo(1:atomNumber),'k','LineWidth',1.5)
plot(meanR1WithLipariSzabo(1:atomNumber),'--k' ...
    ,'LineWidth',1.5)
plot(r1WithSchroedingerEquation(1:atomNumber),'r' ...
    ,'LineWidth',1.5)
plot(meanR1WithSchroedingerEquation(1:atomNumber) ...
    ,'--r','LineWidth',1.5)
legend('Spectral density','Mean Spectral Density' ...
    ,'Lipari Szabo','Mean Lipari Szabo' ...
    ,'Schroedinger Equation','Mean Schroedinger Equ.');
title('Relaxation Rate R1 PSM')
xlabel('Epoches')
ylabel('R1 [Hz]')
grid on
drawnow
hold off
