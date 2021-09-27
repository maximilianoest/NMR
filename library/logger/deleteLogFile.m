function deleteLogFile(path2LogFile)

if isfile(path2LogFile)
    delete(path2LogFile);
end

end
