function plotOrientationDependentResults (path2SaveFigures ...
    ,name2SaveFigures,effectiveMeanRelaxationRates ...
    ,effectiveMedianRelaxationRates,orientations,yLabel)

plottingConfiguration = readConfigurationFile( ...
    'configurationForPlotting.conf');
colors = getColorsFromPlottingConfiguration(plottingConfiguration.colors);
fontSizeOffset = plottingConfiguration.fontSizeOffset;
legendFontSize = plottingConfiguration.legendFontSize+fontSizeOffset;
legendEntries = {};
axisFontSize = plottingConfiguration.axisFontSize+fontSizeOffset;
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
plot(orientations,effectiveMeanRelaxationRates,'Color',colors(1,:))
legendEntries{end+1} = "Mean";
plot(orientations,effectiveMedianRelaxationRates,'Color',colors(2,:))
legendEntries{end+1} = "Median";
lowerBorder = min([min(effectiveMeanRelaxationRates) ...
    min(effectiveMedianRelaxationRates)])*0.95;
upperBorder = max([max(effectiveMeanRelaxationRates) ...
    max(effectiveMedianRelaxationRates)])*1.05;
axis([0 orientations(end) lowerBorder upperBorder]);

legend(legendEntries,'FontSize',legendFontSize,'Interpreter','latex' ...
    ,'Location','east');
legend boxoff 
xlabel('Orientation $\theta \left[ ^{\circ} \right]$')
xticks(0:10:90);
xticklabels(0:10:90);
ylabel(yLabel)

print(gcf,[path2SaveFigures name2SaveFigures '.png'] ,'-dpng'...
    ,resolution)
close all

end
