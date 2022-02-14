function [path2Data,path2Save,path2ConstantsFile,path2LogFile] = ...
    setUpDependenciesBasedOnConfiguration(configuration)

fileName = configuration.fileName;
lipidName = getLipidNameFromFileName(fileName);
constituent = getConstituentFromFileName(fileName);
combinedName = [lipidName constituent];

runOnServer = configuration.runOnServer;
if runOnServer
    path2BaseConfiguration = configuration.path2BaseConfigurationOnServer;
else
    path2BaseConfiguration = ...
        configuration.path2BaseConfigurationOnLocalMachine;
end
    
baseConfiguration = readConfigurationFile(path2BaseConfiguration);
startingDate = datestr(date,'yyyymmdd');

if runOnServer
    path2Data = baseConfiguration.path2DataOnServer;
    path2ConstantsFile = baseConfiguration.path2ConstantsFileOnServer;
    mainPath2Results = baseConfiguration.path2ResultsOnServer;
    additionalFolder = '';
else
    path2Data = baseConfiguration.path2DataOnLocalMachine;
    path2ConstantsFile = baseConfiguration.path2ConstantsFileOnLocalMachine;
    mainPath2Results = baseConfiguration.path2ResultsOnLocalMachine;
    additionalFolder = [configuration.kindOfResults '\'];
end

path2Results = sprintf('%s%s',mainPath2Results,additionalFolder);
if ~isfolder(path2Results)
    mkdir(path2Results);
end

path2Save = sprintf('%s%s_Results_%s_%s.mat',path2Results,startingDate ...
    ,configuration.kindOfResults,combinedName);
path2LogFile = sprintf('%s%s_LogFile_%s%s.txt',path2Results ...
    ,startingDate,configuration.kindOfResults,combinedName); 

end
