function plotOrientationAndPositionDependentResults(path2SaveFigures ...
    ,name2SaveFigures,meanRelaxationRates,positions,orientations,yLabel)

plottingConfiguration = readConfigurationFile( ...
    'configurationForPlotting.conf');
colors = getColorsFromPlottingConfiguration(plottingConfiguration.colors);
fontSizeOffset = plottingConfiguration.fontSizeOffset;
legendFontSize = plottingConfiguration.legendFontSize+fontSizeOffset;
legendEntries = {};
axisFontSize = plottingConfiguration.axisFontSize+fontSizeOffset;
defaultLineWidth = plottingConfiguration.defaultLineWidth;

figure(1)
set(0,'defaultTextInterpreter','latex','DefaultFigureVisible','on');
set(gca,'TickLabelInterpreter','latex')
axisSetUp = gca;
axisSetUp.FontSize = axisFontSize; 
set(0,'DefaultLineLineWidth',defaultLineWidth)
hold on
grid minor
for orientationNr = 1:size(orientations,2)
    plot(positions,meanRelaxationRates(orientationNr,:),'Color' ...
        ,colors(orientationNr,:))
    legendEntries{end+1} = strjoin(["Orientation $\theta$: " ...
        num2str(orientations(orientationNr)) "$^{\circ}$"],''); %#ok<AGROW>
end
lowerBorder = min(min(meanRelaxationRates))*0.95;
upperBorder = max(max(meanRelaxationRates))*1.05;
axis([0 positions(end) lowerBorder upperBorder]);

legend(legendEntries,'FontSize',legendFontSize,'Interpreter','latex');
legend boxoff 
xlabel('Position $\varphi \left[ ^{\circ} \right]$')
xticks(0:10:90);
xticklabels(0:10:90);
ylabel(yLabel)

print(gcf,[path2SaveFigures name2SaveFigures '.png'] ,'-dpng'...
    ,plottingConfiguration.resolution)
close all


end