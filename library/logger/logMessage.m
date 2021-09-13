function logMessage(logInformation,path2LogFile)

fileId = fopen(path2LogFile,'a');
if fileId == -1
    error('Cannot open log file.')
end
fprintf(fileId, '%s: %s\n', datestr(now,0),logInformation);
fclose(fileId);
end
