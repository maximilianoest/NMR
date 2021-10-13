function colors = getStandardizedFigure()

plottingConfiguration = readConfigurationFile( ...
    'configurationForPlotting.txt');
colors = getColorsFromPlottingConfiguration(plottingConfiguration.colors);
fontSizeOffset = plottingConfiguration.fontSizeOffset;

legendFontSize = plottingConfiguration.legendFontSize+fontSizeOffset;
lgd = legend();
lgd.Interpreter = 'latex';
lgd.FontSize = legendFontSize;
legend boxoff

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


end

