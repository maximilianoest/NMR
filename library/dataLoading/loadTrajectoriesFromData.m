function [trajectoryX,trajectoryY,trajectoryZ] = ...
    loadTrajectoriesFromData(configuration,path2Data)

loaded = configuration.dataLoaded;
if not(loaded)
    
    path2File = createMatFilePath(path2Data,configuration.fileName);
    data = load(path2File);
    hydrogenTrajectories = data.(configuration.dataFieldName);
    
    
    trajectoryX = squeeze(hydrogenTrajectories(:,1,:));
    trajectoryY = squeeze(hydrogenTrajectories(:,2,:));
    trajectoryZ = squeeze(hydrogenTrajectories(:,3,:));
else
    
end

end
