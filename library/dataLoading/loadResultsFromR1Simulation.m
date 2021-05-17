function [data] = loadResultsFromR1Simulation(configuration)

compartment = configuration.compartment;
switch compartment
    case "Water"
        fileName = configuration.waterFileName;
    case "Lipid"
        fileName = configuration.lipidFileName;
end
        
data2Load = [configuration.path2Data fileName '.mat'];
data = load(data2Load);

end
