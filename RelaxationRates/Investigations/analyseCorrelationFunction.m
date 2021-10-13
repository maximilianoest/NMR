clc
clear all
close all

% configuration = data.configuration;
configuration = readConfigurationFile('..\config.txt');
addpath(configuration.path2LibraryOnLocalMachine);

path2Data = "C:\Users\maxoe\Google Drive\Promotion\Results\performanceAnalysing\Server\20211008_Results_performanceAnalysing_Lipid_H_500ns_4ps_wh.mat";
data = load(path2Data);
dataConfiguration = data.configuration;
path2Conclusions = [configuration.path2ResultsOnLocalMachine ...
    dataConfiguration.kindOfResults '\Conclusions\'];

correlationFunction0W0 = data.correlationFunction0W0Saver;
correlationFunction1W0 = data.correlationFunction1W0Saver;
correlationFunction2W0 = data.correlationFunction2W0Saver;

normCorrelationFunction0W0 = correlationFunction0W0 ...
    ./correlationFunction0W0(:,:,1);
normCorrelationFunction1W0 = correlationFunction1W0 ...
    ./correlationFunction1W0(:,:,1);
normCorrelationFunction2W0 = correlationFunction2W0 ...
    ./correlationFunction2W0(:,:,1);

dt = data.deltaT;

orientationAngles = rad2deg(data.orientationAngles);
positionAngles = rad2deg(data.positionAngles);
dateNow = datestr(now,'yyyymmdd_HHMM_');

plotCorrelationFunctionsForDifferentOrAndPos(normCorrelationFunction0W0 ...
    ,0,dt,orientationAngles,positionAngles);
saveas(gcf,[path2Conclusions dateNow ...
    'correlationFunctions0W0.png'])
close all

plotCorrelationFunctionsForDifferentOrAndPos(normCorrelationFunction1W0 ...
    ,1,dt,orientationAngles,positionAngles);
saveas(gcf,[path2Conclusions dateNow ...
    'correlationFunctions1W0.png'])
close all

plotCorrelationFunctionsForDifferentOrAndPos(normCorrelationFunction2W0 ...
    ,2,dt,orientationAngles,positionAngles);
saveas(gcf,[path2Conclusions dateNow ...
    'correlationFunctions2W0.png'])
close all

save([path2Conclusions dateNow 'correlationFunctions_Data.mat'])
