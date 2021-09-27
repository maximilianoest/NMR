function [path2Data,path2Save,path2ConstantsFile,path2LogFile] = ...
    setUpSystemBasedOnMachine(configuration)

runOnServer = configuration.runOnServer;
fileName = configuration.fileName;
startingDate = datestr(date,'yyyymmdd');

if runOnServer
    path2Data = configuration.path2DataOnServer;
    path2ConstantsFile = configuration.path2ConstantsFileOnServer;
    path2Results = [configuration.path2ResultsOnServer ...
        configuration.kindOfResults '/'];
else
    path2Data = configuration.path2DataOnLocalMachine;
    path2ConstantsFile = configuration.path2ConstantsFileOnLocalMachine;
    path2Results = [configuration.path2ResultsOnLocalMachine ...
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
