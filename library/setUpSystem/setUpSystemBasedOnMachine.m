function [path2Data,path2Save,path2ConstantsFile] = ...
    setUpSystemBasedOnMachine(configuration)

runOnServer = configuration.runOnServer;

if runOnServer
    path2Data = configuration.path2DataOnServer;
    path2ConstantsFile = configuration.path2ConstantsFileOnServer;
    path2ResultsDirectory = configuration.path2ResultsOnServer; 
else
    path2Data = configuration.path2DataOnLocalMachine;
    path2ConstantsFile = configuration.path2ConstantsFileOnLocalMachine;
    path2ResultsDirectory = configuration.path2ResultsOnLocalMachine;
end

path2Save = [path2ResultsDirectory datestr(date,'yyyymmdd_') fileName ...
    configuration.resultsSuffix];

end
