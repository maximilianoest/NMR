function logMemoryUsage(path2LogFile)
    logMessage(sprintf('Used memory %.6f GB', getUsedMemory),path2LogFile);
end
