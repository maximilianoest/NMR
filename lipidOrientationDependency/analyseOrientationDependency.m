clc; clear all; close all;

results = load("C:\Users\maxoe\Google Drive\Promotion\Results\LIPIDS\DOPS\orientationDependency\20220307_Results_orientationDependency_DOPSlipid.mat");
savingPath = sprintf(['C:\\Users\\maxoe\\Google Drive\\Promotion' ...
    '\\Zwischenergebnisse\\%s_orientationDependency_DOPS\\'] ...
    ,results.startDateOfSimulation);
saving = 1;
if ~isfolder(savingPath)
    mkdir(savingPath);
end
figurePosAndSize = [50 50 900 600];

atomCounter = results.atomCounter;
allR1 = results.scaledR1Rates(:,:,1:atomCounter);
averagedR1 = mean(allR1,3);
orientationAngles = rad2deg(results.orientationAngles);
positionAngles = rad2deg(results.positionAngles);

figure('Position',figurePosAndSize);
legendEntries = {};
set(gcf,'DefaultLineLineWidth',1.5)
set(gca,'FontSize',16)
hold on
for positionNr = 1:length(positionAngles)
    legendEntries{end+1} = sprintf('$\\varphi$ = %.1f$^{\\circ}$' ...
        ,positionAngles(positionNr));
    p = plot(orientationAngles ...
        ,squeeze(averagedR1(:,positionNr)));
end
hold off
grid minor
legend(legendEntries,'Interpreter','latex','Location','northwest');
title('DOPS orientation-dependent R1','Interpreter','latex');
xlabel('Orientation angle $\theta$ [$^{\circ}$]','Interpreter','latex');
ylabel('Relaxation rate [Hz]','Interpreter','latex');
if saving
    savingName = sprintf('orientationDependentR1_%s' ...
        ,results.whichLipid); %#ok<UNRCH>
    print(gcf,[savingPath savingName],'-dpng','-r300');
end

effectiveR1 = squeeze(mean(averagedR1,2));
figure('Position',figurePosAndSize);
legendEntries = {};
set(gcf,'DefaultLineLineWidth',1.5)
set(gca,'FontSize',16)
hold on
p = plot(orientationAngles,effectiveR1);
hold off
grid minor
title('DOPS orientation-dependent effective R1','Interpreter','latex');
xlabel('Orientation angle $\theta$ [$^{\circ}$]','Interpreter','latex');
ylabel('Effective relaxation rate [Hz]','Interpreter','latex');
if saving
    savingName = sprintf('orientationDependentEffectiveR1_%s' ...
        ,results.whichLipid); %#ok<UNRCH>
    print(gcf,[savingPath savingName],'-dpng','-r300');
end


