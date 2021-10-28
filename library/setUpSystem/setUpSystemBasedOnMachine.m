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

fileNames = getValuesFromStringEnumeration(configuration.fileNames,';' ...
    ,'string');
fileNameArray = string();
for fileNameNr = 1:length(fileNames)
    fileName = fileNames(fileNameNr);
    fileName = strsplit(fileName,'_');
    fileName = fileName(1);
    if ~contains(fileNameArray,fileName)
        fileNameArray = fileNameArray + fileName;
    end
end
fileNameArray = convertStringsToChars(fileNameArray);

path2Save = [path2Results startingDate '_Results_' ...
    configuration.kindOfResults '_' fileNameArray '.mat'];
path2LogFile = [path2Results startingDate '_LogFile_' ...
    configuration.kindOfResults '_' fileNameArray '.txt'];

end
