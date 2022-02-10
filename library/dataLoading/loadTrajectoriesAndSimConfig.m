function [trajectoryX,trajectoryY,trajectoryZ,simConfig] ...
    = loadTrajectoriesAndSimConfig(path2File)

matObject = matfile(path2File);
simConfig = matObject.configuration;

load(path2File,'trajectories');

trajectoryX = squeeze(single(trajectories(:,1,:)));
trajectoryY = squeeze(single(trajectories(:,2,:)));
trajectoryZ = squeeze(single(trajectories(:,3,:)));

end