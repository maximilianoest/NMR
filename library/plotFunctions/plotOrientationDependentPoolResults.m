function plotOrientationDependentPoolResults(path2SaveFigures ...
    ,name2SaveFigures,predictedR1Rates,orientations,yLabel)

plottingConfiguration = readConfigurationFile( ...
    'configurationForPlotting.conf');
colors = getColorsFromPlottingConfiguration(plottingConfiguration.colors);
legendFontSize = plottingConfiguration.legendFontSize;
legendEntries = {};
axisFontSize = plottingConfiguration.axisFontSize;
defaultLineWidth = plottingConfiguration.defaultLineWidth;
resolution = plottingConfiguration.resolution;

figure(1)
set(0,'defaultTextInterpreter','latex','DefaultFigureVisible','on');
set(gca,'TickLabelInterpreter','latex')
axisSetUp = gca;
axisSetUp.FontSize = axisFontSize; 
set(0,'DefaultLineLineWidth',defaultLineWidth)
hold on
grid minor

plot(orientations,predictedR1Rates(1,:),':','Color',colors(1,:))
legendEntries{end+1} = "Case 1";
plot(orientations,predictedR1Rates(2,:),'--','Color',colors(1,:))
legendEntries{end+1} = "Case 2";
plot(orientations,predictedR1Rates(3,:),':','Color',colors(2,:))
legendEntries{end+1} = "Case 3";
plot(orientations,predictedR1Rates(4,:),'--','Color',colors(2,:))
legendEntries{end+1} = "Case 4";
lowerBorder = min(min(predictedR1Rates))*0.999;
upperBorder = max(max(predictedR1Rates))*1.001;
axis([0 orientations(end) lowerBorder upperBorder]);

legend(legendEntries,'FontSize',legendFontSize,'Interpreter','latex' ...
    ,'Location','northwest');
legend boxoff 
xlabel('Orientation $\theta \left[ ^{\circ} \right]$')
xticks(0:10:90);
xticklabels(0:10:90);
ylabel(yLabel)

print(gcf,[path2SaveFigures name2SaveFigures '.png'] ,'-dpng'...
    ,resolution)
close all

end