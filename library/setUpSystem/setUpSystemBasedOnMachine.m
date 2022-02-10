function [path2Data,path2Save,path2ConstantsFile,path2LogFile] = ...
    setUpSystemBasedOnMachine(configuration)

error('This function is old, please use setUpSystemBasedOnConfiguration');
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

if runOnServer
    addtionalFolder = '';
else   
    additionalFolder = [configuration.kindOfResults '_'];
    
end
path2Save = sprintf('%s%s_Results_%s%s.mat',path2Results,startingDate ...
    ,additionalFolder,fileName);
path2LogFile = sprintf('%s%s_LogFile_%s%s.txt',path2Results ...
    ,startingDate,additionalFolder,fileName); 

end
