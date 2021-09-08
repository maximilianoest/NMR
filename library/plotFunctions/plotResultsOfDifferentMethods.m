function plotResultsOfDifferentMethods(path2SaveFigures ...
    ,name2SaveFigures,path2Data,folderFileName,numberOfEpochs)

data = load([path2Data folderFileName '.mat']);
meanResultsFromMethods = [data.meanR1WithPerturbationTheory' ...
    ,data.meanR1WithLipariSzabo',data.meanR1WithSchroedingerEquation']';

r1DataFromMethods = [data.r1WithPerturbationTheory' ...
    ,data.r1WithLipariSzabo',data.r1WithSchroedingerEquation']';

plottingConfiguration = readConfigurationFile( ...
    'configurationForPlotting.conf');
colors = getColorsFromPlottingConfiguration(plottingConfiguration.colors);
fontSizeOffset = plottingConfiguration.fontSizeOffset;
legendFontSize = plottingConfiguration.legendFontSize + fontSizeOffset;
axisFontSize = plottingConfiguration.axisFontSize + fontSizeOffset;
defaultLineWidth = plottingConfiguration.defaultLineWidth;
resolution = plottingConfiguration.resolution;

figure(1)
fig = gcf;
fig.PaperPositionMode = 'auto';
set(0,'defaultTextInterpreter','latex','DefaultFigureVisible','on');
set(gca,'TickLabelInterpreter','latex')
axisSetUp = gca;
axisSetUp.FontSize = axisFontSize; 
set(0,'DefaultLineLineWidth',defaultLineWidth)
grid minor
hold on
disp(name2SaveFigures)
disp('All three relaxation rates at the end: ')
methods = ["Explicit FT: " "Lipari Szabo: " "Schroedinger: "];
for methodNr = 1:3
    results = meanResultsFromMethods(methodNr,1:numberOfEpochs);
    plot(results,'Color',colors(methodNr,1:3));
    disp(strjoin(["    " methods(methodNr) num2str(results(end)) ...
        " \pm " num2str(std(r1DataFromMethods(methodNr,:))) " Median: " ...
        num2str(median(results))]))
end
upperBorder = ceil(max(max(meanResultsFromMethods(1:numberOfEpochs))))+1;
lowerBorder = floor(min(min(meanResultsFromMethods(1:numberOfEpochs))))-1;
axis([0 numberOfEpochs lowerBorder upperBorder]);
legend({'Explicit FT','Lipari-Szabo','Schr\"odinger'} ...
    ,'FontSize',legendFontSize,'Interpreter','latex');
legend boxoff 
xlabel('Number of Hydrogen Nuclei')
xticks(0:20:140);
xticklabels(0:20:140);
ylabel('Relaxation Rate $R_1 \left[ Hz \right]$')

print(gcf,[path2SaveFigures name2SaveFigures '.png'] ,'-dpng'...
    ,resolution)
hold off
close all


end