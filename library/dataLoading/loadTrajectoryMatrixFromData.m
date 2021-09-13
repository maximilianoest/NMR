function [hydrogenTrajectories] = ...
    loadTrajectoryMatrixFromData(configuration,path2Data)

loaded = configuration.dataLoaded;
if not(loaded)
    
    fileName = configuration.fileName;
    path2File = [path2Data fileName];
    data = load(path2File);
    hydrogenTrajectories = data.(configuration.dataFieldName);
    
else
    
end

end
