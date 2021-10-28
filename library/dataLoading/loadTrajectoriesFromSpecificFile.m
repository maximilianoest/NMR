function [trajectoryX,trajectoryY,trajectoryZ] = ...
    loadTrajectoriesFromSpecificFile(configuration,path2Data,fileName ...
    ,fieldName)

loaded = configuration.dataLoaded;
if not(loaded)
    
    path2File = createMatFilePath(path2Data,fileName);
    data = load(path2File);
    hydrogenTrajectories = data.(fieldName);
     
    trajectoryX = single(squeeze(hydrogenTrajectories(:,1,:)));
    trajectoryY = single(squeeze(hydrogenTrajectories(:,2,:)));
    trajectoryZ = single(squeeze(hydrogenTrajectories(:,3,:)));
else
    
end

end
