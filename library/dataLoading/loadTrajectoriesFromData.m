function [trajectoryX,trajectoryY,trajectoryZ] = ...
    loadTrajectoriesFromData(configuration,path2Data)

loaded = configuration.dataLoaded;
if not(loaded)
    
    path2File = createMatFilePath(path2Data,configuration.fileName);
    data = load(path2File);
    hydrogenTrajectories = data.(configuration.dataFieldName);
    
    
    trajectoryX = single(squeeze(hydrogenTrajectories(:,1,:)));
    trajectoryY = single(squeeze(hydrogenTrajectories(:,2,:)));
    trajectoryZ = single(squeeze(hydrogenTrajectories(:,3,:)));
else
    
end

end
