function printBreakLineToLogFile(path2LogFile)

fileId = fopen(path2LogFile,'a');
if fileId == -1
    error('Cannot open log file.')
end
fprintf(fileId,'____________________________________________________\n');
fclose(fileId);
end
