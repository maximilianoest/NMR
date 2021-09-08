function plotCorrelationFunctionsOfConstituents(path2SaveFigures ...
    ,name2SaveFigures,path2Data,folderFileName,dt,colorNr)

data = load([path2Data folderFileName '.mat']);
normalizedCorrelationFunctions = data.normCorrelation;

plottingConfiguration = readConfigurationFile( ...
    'configurationForPlotting.conf');
colors = getColorsFromPlottingConfiguration(plottingConfiguration.colors);
timeAxis = 0:dt:size(normalizedCorrelationFunctions,2)*dt-dt;

disp(name2SaveFigures)

plot(timeAxis,normalizedCorrelationFunctions(1,:) ...
    ,'--','Color',colors(colorNr,:))
plot(timeAxis,normalizedCorrelationFunctions(2,:) ...
    ,':','Color',colors(colorNr,:))

upperBorderY = 0.4;
lowerBorderY = 0;
upperBorderX = timeAxis(end);
lowerBorderX = timeAxis(1);
axis([lowerBorderX upperBorderX lowerBorderY upperBorderY]);


end