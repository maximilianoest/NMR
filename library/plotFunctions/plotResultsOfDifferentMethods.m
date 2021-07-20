function plotResultsOfDifferentMethods(path2SaveFigures ...
    ,name2SaveFigures,path2Data,folderFileName,numberOfEpochs ...
    ,resolution)

data = load([path2Data folderFileName '.mat']);
meanResultsFromMethods = [data.meanR1WithPerturbationTheory' ...
    ,data.meanR1WithLipariSzabo',data.meanR1WithSchroedingerEquation']';

% Convert color code to 1-by-3 RGB array (0~1 each)
explicitFTColor = '#6699ff';
explicitFTColor = sscanf(explicitFTColor(2:end),'%2x%2x%2x',[1 3])/255;
liSzColor = '#00cc66';
liSzColor = sscanf(liSzColor(2:end),'%2x%2x%2x',[1 3])/255;
schroedingerColor = '#ff9933';
schroedingerColor = sscanf(schroedingerColor(2:end),'%2x%2x%2x',[1 3])/255;
colors = [explicitFTColor',liSzColor',schroedingerColor']';

figure(1)
legendFontSize = 15;
fig = gcf;
fig.PaperPositionMode = 'auto';
set(0,'defaultTextInterpreter','latex','DefaultFigureVisible','on');
set(gca,'TickLabelInterpreter','latex')
axisSetUp = gca;
axisSetUp.FontSize = 15; 
set(0,'DefaultLineLineWidth',1.7)
grid minor
hold on
disp(name2SaveFigures)
disp('All three relaxation rates at the end: ')
methods = ["Explicit FT: " "Lipari Szabo: " "Schroedinger: "];
for methodNr = 1:3
    results = meanResultsFromMethods(methodNr,1:numberOfEpochs);
    plot(results,'Color',colors(methodNr,1:3));
    disp(strjoin(["    " methods(methodNr) num2str(results(end)) ...
        " \pm " num2str(std(results)) " Median: " ...
        num2str(median(results))]))
end
upperBorder = ceil(max(max(meanResultsFromMethods(1:numberOfEpochs))))+1;
lowerBorder = floor(min(min(meanResultsFromMethods(1:numberOfEpochs))))-1;
axis([0 numberOfEpochs lowerBorder upperBorder]);
legend({'Explicit FT','Lipari-Szabo','Schr\"odinger'} ...
    ,'FontSize',legendFontSize,'Interpreter','latex');
legend boxoff 
xlabel('Number of Epochs')
ylabel('Relaxation Rate [Hz]')

print(gcf,[path2SaveFigures name2SaveFigures '.png'] ,'-dpng'...
    ,resolution)
hold off
close all


end