function [trajectoryX,trajectoryY,trajectoryZ] = ...
    loadTrajectoriesFromData(configuration,path2Data)

loaded = configuration.dataLoaded;
if not(loaded)
    
    disp('Loading data')
    
    fileName = configuration.fileName;
    path2File = [path2Data fileName '.mat'];
    data = load(path2File);
    hydrogenTrajectories = data.(configuration.dataFieldName);
    
    disp('Data successfully loaded')
    
    trajectoryX = squeeze(hydrogenTrajectories(:,1,:));
    trajectoryY = squeeze(hydrogenTrajectories(:,2,:));
    trajectoryZ = squeeze(hydrogenTrajectories(:,3,:));
    
end

end
