clc; clear all; close all;

results = load("C:\Users\maxoe\Google Drive\Promotion\Results\LIPIDS\PLPC\necessaryNearestNeighbours\20220419_Results_relevantNearestNeighbours_PLPClipid.mat");

r1 = results.r1WithPerturbationTheory;
saving = 1;
savingPath = sprintf(['C:\\Users\\maxoe\\Google Drive\\Promotion' ...
    '\\Zwischenergebnisse\\%s_numberOfAtoms_%s\\'] ...
    ,datestr(date,'yyyymmdd'),results.whichLipid);
if ~exist(savingPath,'dir')
    mkdir(savingPath);
end


orientationAngles = results.orientationAngles;
positionAngles = results.positionAngles;

atomNrStepsSize = 10;
atomCounter = round(results.atomCounter/atomNrStepsSize);
averageR1 = zeros(1,atomCounter);
stdR1 = zeros(1,atomCounter);
numberOfAtomsAxis = 1:atomNrStepsSize:results.atomCounter;

figurePosAndSize = [50 50 900 600];
figure('Position',figurePosAndSize);
set(gcf,'DefaultLineLineWidth',1.5)
set(gca,'FontSize',16)

legendEntries = {};
hold on
for orientationNr = 1:length(orientationAngles)
    for positionNr = 1:length(positionAngles)
       r1ForOriAndPos = squeeze(r1(orientationNr,positionNr,:,1));
        for atomNr = 1:atomCounter
           averageR1(atomNr) = mean( ...
               r1ForOriAndPos(1:(atomNr-1)*atomNrStepsSize+1));
           stdR1(atomNr) = std(...
               r1ForOriAndPos(1:(atomNr-1)*atomNrStepsSize+1));
        end 
       errorbar(numberOfAtomsAxis,averageR1,stdR1,'LineWidth', 1.3);
       legendEntries{end+1} = sprintf("$\\theta$: %.2f, $\\varphi$: %.2f" ...
           ,rad2deg(orientationAngles(orientationNr)) ...
           ,rad2deg(positionAngles(positionNr))); %#ok<SAGROW>
    end
end
hold off
legend(legendEntries,'Interpreter','latex');
xlabel('Number of calculated atoms');
ylabel('Relaxation Rate with STD [Hz]');
title(sprintf('Relaxation rate in dependence of calculated atoms (%s)' ...
    ,results.whichLipid))
grid minor

if saving
    savingName = sprintf('r1DependentOnCalculatedAtoms_%s' ...
        ,results.whichLipid); %#ok<UNRCH>
    print(gcf,[savingPath savingName],'-dpng','-r300');
end
