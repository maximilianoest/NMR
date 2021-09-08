function plotIntraVsInterDependentResults(waterIntraMeanShift ...
    ,waterIntraMedianShift,waterInterMeanShift,waterInterMedianShift ...
    ,path2SaveFigures,name2SaveFigures,yLabel,orientations)

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
plot(orientations,waterIntraMeanShift,':','Color' ...
    ,colors(1,:))
legendEntries{end+1} = "Intra Mean";
plot(orientations,waterIntraMedianShift,'--','Color' ...
    ,colors(1,:))
legendEntries{end+1} = "Intra Median";
plot(orientations,waterInterMeanShift,':','Color' ...
    ,colors(2,:))
legendEntries{end+1} = "Intra + Extra Mean";
plot(orientations,waterInterMedianShift,'--','Color' ...
    ,colors(2,:))
legendEntries{end+1} = "Intra + Extra Median";

hold off

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