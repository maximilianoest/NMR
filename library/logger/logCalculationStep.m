function logCalculationStep(message,path2LogFile)

fileId = fopen(path2LogFile,'a');
if fileId == -1
    error('Cannot open log file.')
end
fprintf(fileId,'%s: %.6f GB \n', message,getUsedMemory());
fclose(fileId);

end
