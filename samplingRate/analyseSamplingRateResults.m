clc
clear all
close all
%#ok<*SAGROW>

results = load("C:\Users\maxoe\Google Drive\Promotion\Results\LIPIDS\DOPS\samplingRates\20220228_Results_samplingRateOptimization_DOPSlipid.mat");
savingPath = ['C:\Users\maxoe\Google Drive\Promotion' ...
    '\Zwischenergebnisse\20220228_samplingRates_DOPS\'];
saving = 1;

figurePosAndSize = [50 50 900 600];

allR1 = results.scaledR1Rates(:,:,1:results.atomCounter,:);
averageR1 = squeeze(mean(allR1,3));
samplingRates = results.samplingFrequencyArray;
calculatedAtoms = results.atomsToCalculate;

orientationAngles = rad2deg(results.orientationAngles);
positionAngles = rad2deg(results.positionAngles);

lineStyle = ["--",":","-."];


for positionNr = 1:length(positionAngles)
    figure('Position',figurePosAndSize);
    legendEntries = {};
    set(gcf,'DefaultLineLineWidth',1.5)
    set(gca,'FontSize',16)
    hold on
    for orientationNr = 1:1:length(orientationAngles)
        legendEntries{end+1} = sprintf('$\\theta$ = %.0f$^{\\circ}$' ...
            ,orientationAngles(orientationNr));
        p = plot(samplingRates./1e-12 ...
            ,squeeze(averageR1(orientationNr,positionNr,:)));
%         p.LineStyle = lineStyle(positionNr);
    end
    hold off
    grid minor
    legend(legendEntries,'Interpreter','latex','Location','northwest')
    title(sprintf(['DOPS sampling rate dependent R1 ' ...
        ', Position $\\varphi$: %.1f$^{\\circ}$'] ...
        ,positionAngles(positionNr)),'Interpreter','latex');
    xlabel('Time between two time steps [ps]','Interpreter','latex')
    ylabel('Relaxation rate [Hz]','Interpreter','latex')
    if saving
        savingName = sprintf('samplingRateDependentR1_%s_Position%i' ...
            ,results.whichLipid,positionAngles(positionNr)); %#ok<UNRCH>
        print(gcf,[savingPath savingName],'-dpng','-r300');
    end
end

effectiveR1 = squeeze(mean(averageR1,2));
figure('Position',figurePosAndSize);
legendEntries = {};
set(gcf,'DefaultLineLineWidth',1.5)
set(gca,'FontSize',16)
hold on

for orientationNr = 1:length(orientationAngles)
    legendEntries{end+1} = sprintf('$\\theta$ = %.0f$^{\\circ}$' ...
        ,orientationAngles(orientationNr));
    p = plot(samplingRates./1e-12 ...
        ,squeeze(effectiveR1(orientationNr,:)));
end

hold off
grid minor
legend(legendEntries,'Interpreter','latex','Location','northwest')
title(sprintf(['DOPS sampling rate dependent effective R1 '] ...
    ,positionAngles(positionNr)),'Interpreter','latex');
xlabel('Time between two time steps [ps]','Interpreter','latex')
ylabel('Effective relaxation rate [Hz]','Interpreter','latex')
if saving
    savingName = sprintf('samplingRateDependentEffectiveR1_%s' ...
        ,results.whichLipid); %#ok<UNRCH>
    print(gcf,[savingPath savingName],'-dpng','-r300');
end

effectiveT1 = 1./squeeze(mean(averageR1,2));
figure('Position',figurePosAndSize);
legendEntries = {};
set(gcf,'DefaultLineLineWidth',1.5)
set(gca,'FontSize',16)
hold on

for orientationNr = 1:length(orientationAngles)
    legendEntries{end+1} = sprintf('$\\theta$ = %.0f$^{\\circ}$' ...
        ,orientationAngles(orientationNr));
    p = plot(samplingRates./1e-12 ...
        ,squeeze(effectiveT1(orientationNr,:)));
end

hold off
grid minor
legend(legendEntries,'Interpreter','latex','Location','northwest')
title(sprintf(['DOPS sampling rate dependent effective T1 '] ...
    ,positionAngles(positionNr)),'Interpreter','latex');
xlabel('Time between two time steps [ps]','Interpreter','latex')
ylabel('Effective relaxation time [sec]','Interpreter','latex')
% if saving
%     savingName = sprintf('samplingRateDependentEffectiveR1_%s' ...
%         ,results.whichLipid); %#ok<UNRCH>
%     print(gcf,[savingPath savingName],'-dpng','-r300');
% end



