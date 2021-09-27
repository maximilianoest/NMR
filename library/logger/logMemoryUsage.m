function logMemoryUsage(path2LogFile)
configuration = readConfigurationFile('config.txt');
if configuration.runOnServer
    usedMemory = getUsedMemoryOnLinux;
else
    usedMemory = getUsedMemoryOnWindows; 
end
logMessage(sprintf('    Used memory %.3f MB', usedMemory) ...
    ,path2LogFile,false);

end
