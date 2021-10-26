function [path2Data,path2Save,path2ConstantsFile,path2LogFile] = ...
    setUpSystemBasedOnMachine(configuration)

runOnServer = configuration.runOnServer;
if runOnServer
    path2BaseConfiguration = configuration.path2BaseConfigurationOnServer;
else
    path2BaseConfiguration = ...
        configuration.path2BaseConfigurationOnLocalMachine;
end
    
baseConfiguration = readConfigurationFile(path2BaseConfiguration);
fileName = configuration.fileName;
startingDate = datestr(date,'yyyymmdd');

if runOnServer
    path2Data = baseConfiguration.path2DataOnServer;
    path2ConstantsFile = baseConfiguration.path2ConstantsFileOnServer;
    path2Results = [baseConfiguration.path2ResultsOnServer ...
        configuration.kindOfResults '/'];
else
    path2Data = baseConfiguration.path2DataOnLocalMachine;
    path2ConstantsFile = baseConfiguration.path2ConstantsFileOnLocalMachine;
    path2Results = [baseConfiguration.path2ResultsOnLocalMachine ...
        configuration.kindOfResults '\'];
end

if ~isfolder(path2Results)
    mkdir(path2Results);
end

path2Save = [path2Results startingDate '_Results_' ...
    configuration.kindOfResults '_' fileName '.mat'];
path2LogFile = [path2Results startingDate '_LogFile_' ...
    configuration.kindOfResults '_' fileName '.txt'];

end
