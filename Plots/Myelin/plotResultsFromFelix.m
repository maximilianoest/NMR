clc 
clear all
close all

data = load(['C:\Users\maxoe\Google Drive\Promotion\Results\FELIX' ...
    '\Results_Lipid_AD_250ns_2ps_FINAL.mat']);

relaxationRateACF = mean(data.r1_ACF,1);
effectiveRelaxationRateACF = mean(data.R1_ACF,1);

relaxationRateSG = mean(data.r1_SG,1);
effectiveRelaxationRateSG = mean(data.R1_SG,1);

figure(1)
plot(effectiveRelaxationRateACF)
hold on
plot(effectiveRelaxationRateSG)
hold off
legend('ACF','SG')
grid on


% r1 = data.R1;
% meanR1 = mean(r1,1);
% theta = rad2deg(data.Theta);
% meanR1 = meanR1(1:13);
% difference = meanR1(end)- meanR1(1)
% 
% figure(1)
% plot(meanR1)


