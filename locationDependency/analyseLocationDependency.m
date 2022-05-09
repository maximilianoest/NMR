clc; clear all; close all


%% initialize
results = load("C:\Users\maxoe\Google Drive\Promotion\Results\LIPIDS\PLPC\orientationDependency\20220502_Results_orientationDependency_PLPClipid.mat");
whichLipid = results.whichLipid;
savingPath = initializeSystemForSavingPlots("orientationDependency" ...
    ,whichLipid);
saving = 1;
trimming = 0;
trimmingString = ''; %#ok<NASGU>

atomCounter = results.atomCounter;
allR1 = results.scaledR1Rates(:,:,1:atomCounter);
averagedR1 = mean(allR1,3);
orientationAngles = rad2deg(results.orientationAngles);
positionAngles = rad2deg(results.positionAngles);
orientationCount = size(orientationAngles,2);
positionsCount = size(positionAngles,2);


%% trimm data
if trimming
    lowerPercentile = 10; %#ok<*UNRCH>
    upperPercentile = 80;

    trimmedR1 = trimmR1DataWithPercentile(allR1,lowerPercentile ...
        ,upperPercentile,orientationCount,positionsCount);
    trimmingString = sprintf('TRIMMED_%.0f_%.0f',lowerPercentile ...
        ,upperPercentile);
    allR1 = trimmedR1;
end

%% determine STD
standardDeviations = zeros(orientationCount,positionsCount)*NaN;
for orientationNr = 1:orientationCount
    for positionNr = 1:positionsCount
        standardDeviations(orientationNr,positionNr) = ...
            std(allR1(orientationNr,positionNr,:));
    end
end

%% plot and save stuff
initializeFigure();
legendEntries = {};
for positionNr = 1:length(positionAngles)
    legendEntries{end+1} = sprintf('$\\varphi$ = %.1f$^{\\circ}$' ...
        ,positionAngles(positionNr)); %#ok<SAGROW>
    plot(orientationAngles ...
        ,squeeze(averagedR1(:,positionNr)));
end
legend(legendEntries,'Location','northwest');
title(sprintf('Orientation-dependent R1 (%s)',whichLipid));
xlabel('Orientation angle $\theta$ [$^{\circ}$]');
ylabel('Relaxation rate [Hz]');
if saving
    savingName = sprintf('%s_%sorientationDependentR1_%s' ...
        ,results.startDateOfSimulation,trimmingString,whichLipid); 
    print(gcf,[savingPath savingName],'-dpng','-r300');
end

initializeFigure('legend',false);
averagedR1 = mean(allR1,3);
orientationAngles = rad2deg(results.orientationAngles);
positionAngles = rad2deg(results.positionAngles);
positionAxis = [];
for orientationNr = 1:size(orientationAngles,2)
    positionAxis = [positionAxis;positionAngles]; %#ok<AGROW>
end
orientationAxis = [];
for positionNr = 1:size(positionAngles,2)
    orientationAxis = [orientationAxis;orientationAngles]; %#ok<AGROW>
end
orientationAxis = orientationAxis';
coloring(:,:,1) = zeros(3,5); % red
coloring(:,:,2) = ones(3,5).*linspace(0.2,0.8,3)'; % green
coloring(:,:,3) = ones(3,5).*linspace(0.2,0.8,5); % blue

surf(positionAxis,orientationAxis,averagedR1,coloring,'FaceAlpha',0.8);
viewAngle = -45;
additionalAngle = 30;
view(viewAngle,30)
title(sprintf('Orientation- and position-dependent R1 (%s)' ...
    ,whichLipid));
xlabel('$\varphi$ [$^{\circ}$]')%,'rotation', -viewAngle-additionalAngle);
ylabel('$\theta$ [$^{\circ}$]')%,'rotation', viewAngle+additionalAngle);
zlabel('Relaxation rate [Hz]');
if saving
    savingName = sprintf('%s_%s3dDependentR1_%s' ...
        ,results.startDateOfSimulation,trimmingString,whichLipid); 
    print(gcf,[savingPath savingName],'-dpng','-r300');
end


effectiveR1 = squeeze(mean(averagedR1,2)) %#ok<NOPTS>
initializeFigure('legend',false);
p = plot(orientationAngles,effectiveR1);
title(sprintf('Orientation-dependent effective R1 (%s)',whichLipid));
xlabel('Orientation angle $\theta$ [$^{\circ}$]');
ylabel('Effective relaxation rate [Hz]');

overallR1 = mean(effectiveR1) %#ok<NOPTS>

if saving
    savingName = sprintf('%s_%sorientationDependentEffectiveR1_%s' ...
        ,results.startDateOfSimulation,trimmingString,results.whichLipid); 
    print(gcf,[savingPath savingName],'-dpng','-r300');
    formatToSave = [orientationAngles; effectiveR1'];
    fileId = fopen([savingPath sprintf('%s_orientationDependentR1_%s.txt' ...
        ,results.startDateOfSimulation,results.whichLipid)],'a');
    if trimming
        fprintf(fileId,['TRIMMED DATA (percentiles lower = %.2f, upper ' ...
            '= %.2f) \r\n'],lowerPercentile,upperPercentile);
    else
        fprintf(fileId,'UNTRIMMED DATA: \r\n');
    end
    fprintf(fileId,'Relaxation Rate (Standard deviation): \r\n');
    fprintf(fileId,'%11s','Theta\\Phi');
    fprintf(fileId,'%20.2f',positionAngles);
    fprintf(fileId,'\n');
    for orientationNr = 1:orientationCount
        fprintf(fileId,'%11.2f',orientationAngles(orientationNr));
        for positionNr = 1:positionsCount
            fprintf(fileId,'   %7.4f (%7.4f)',averagedR1(orientationNr ...
                ,positionNr),standardDeviations(orientationNr,positionNr));
        end
        fprintf(fileId,'\r\n');
    end
    fprintf(fileId,'\r\n\r\n');
    fprintf(fileId,'Effective R1: \r\n');
    fprintf(fileId,'\tTheta = %6.2f: %.4f\r\n',formatToSave);
    fprintf(fileId,'\r\nOverall R1: \r\n');
    fprintf(fileId,'\t%10s%.4f \r\n',' ',overallR1);
    fprintf(fileId,'----------------------------------------- \r\n\r\n');
    fclose(fileId);
end
