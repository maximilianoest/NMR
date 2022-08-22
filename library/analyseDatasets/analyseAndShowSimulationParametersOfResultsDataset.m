function analyseAndShowSimulationParametersOfResultsDataset(dataset)
% Load dataset and present data. 

fprintf('Some information about the dataset:\n');
filePAth = dataset.path2File;
fprintf('  File path: %s\n',filePAth);
numberOfCalculatedAtoms = dataset.atomCounter;
fprintf('  Gromacs simulation temperature: %.2f K\n', ...
    table(str2num(dataset.simulationConfiguration.ref_t)).Var1(1));
fprintf('  Number of calculated atoms: %i\n',numberOfCalculatedAtoms);
numberOfNearestNeighbours = dataset.nearestNeighbours;
fprintf('  Number of nearest neighbours: %i \n', ...
    numberOfNearestNeighbours);
thetaAngles = rad2deg(dataset.orientationAngles);
fprintf(['  Theta angles:'  repmat(' %.2f',1,length(thetaAngles)), ...
    '\n'],thetaAngles);
positionAngles = rad2deg(dataset.positionAngles);
fprintf(['  Phi angles:'  repmat(' %.2f',1,length(positionAngles)), ...
    '\n'],positionAngles);


end