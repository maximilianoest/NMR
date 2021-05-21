%% Set Up System
clc

configuration = readConfigurationFile('config.conf');

if configuration.runOnServer
    addpath(genpath(configuration.path2LibraryOnServer));
else
    addpath(genpath(configuration.path2LibraryOnLocalMachine))
end

calculateR1RateRAMandPERFORMANCEoptimized(configuration);


