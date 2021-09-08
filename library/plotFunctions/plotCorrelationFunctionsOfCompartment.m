function plotCorrelationFunctionsOfCompartment(path2SaveFigures ...
    ,name2SaveFigures,correlationFunction1W0,correlationFunction2W0 ...
    ,yLabel,timeAxis)

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
set(gca,'XMinorTick','on','YMinorTick','on')
axisSetUp = gca;
axisSetUp.FontSize = axisFontSize; 
set(0,'DefaultLineLineWidth',defaultLineWidth)
hold on
grid minor

plot(timeAxis,correlationFunction1W0,'Color' ...
    ,colors(1,:))
legendEntries{end+1} = '$q = 1$';
plot(timeAxis,correlationFunction2W0,'Color' ...
    ,colors(2,:))
legendEntries{end+1} = '$q = 2$';

axis([0 timeAxis(end) 0 0.5]);

legend(legendEntries,'FontSize',legendFontSize,'Interpreter','latex');
legend boxoff 
xlabel('$\tau [s]$')
ylabel(yLabel)
% yticks(0:0.1:1)
% yticklabels(0:0.1:1)

print(gcf,[path2SaveFigures name2SaveFigures '.png'] ,'-dpng'...
    ,plottingConfiguration.resolution)
close all

end