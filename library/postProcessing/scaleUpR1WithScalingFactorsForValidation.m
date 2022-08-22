function scaledUpR1_path = scaleUpR1WithScalingFactorsForValidation( ...
    r1ForDifferentNN_path,scalFac_path)

% scaling factors
scalFacResults = load(scalFac_path);
scalFac_NN = scalFacResults.scalFac_NN;


% relaxation rates to scale up
r1ForDifferentNN = load(r1ForDifferentNN_path);
r1_theta_phi_NN = r1ForDifferentNN.r1_theta_phi_NN;

% check for compatibility
testScalingFactorsAndR1sToScaleForCompatibility(scalFacResults, ...
    r1ForDifferentNN);

% define parameters
whichLipid = r1ForDifferentNN.whichLipid;
fprintf('Lipid: %s\n',whichLipid);
fieldStrength = r1ForDifferentNN.fieldStrength;
fprintf('Field strength: %.3f\n',fieldStrength);
scalFacMatlabSimulationDate = scalFacResults.matlabSimulationDate;
fprintf('Start date simulation scal fac: %s\n',scalFacMatlabSimulationDate);
r1MatlabSimulationDate = r1ForDifferentNN.matlabSimulationDate;
fprintf('Start date simulation R1: %s\n',r1MatlabSimulationDate);
scalFacGromacsSimulationDate = scalFacResults.gromacsSimulationDate;
r1GromascSimulationDate = r1ForDifferentNN.gromacsSimulationDate;


fprintf('NN cases scal Fac: ');
fprintf('%i  ',scalFacResults.nearestNeighbourCases);
fprintf('\nNN cases R1:     ');
nearestNeighbourCases = r1ForDifferentNN.nearestNeighbourCases;
fprintf('%i  ',nearestNeighbourCases);
fprintf('\n');
orientationAngles = r1ForDifferentNN.orientationAngles;
positionAngles = r1ForDifferentNN.positionAngles;

% scaling up with scaling factor SF(NN)
for oriNr = 1:length(orientationAngles)
    for posNr = 1:length(positionAngles)
        r1ScaledUp_theta_phi_NN(oriNr,posNr,:) = ...
            squeeze(r1_theta_phi_NN(oriNr,posNr,:))'.*scalFac_NN; %#ok<AGROW>
    end
end
r1ScaledUp_theta_NN = squeeze(mean(r1ScaledUp_theta_phi_NN,2));
r1ScaledUp_NN = squeeze(mean(r1ScaledUp_theta_NN,1)); %#ok<NASGU>

savingPath = initializeSystemForSavingR1();
fieldstrengthString = strrep(num2str(fieldStrength),'.','');
savingName = sprintf('%s_r1%s_scalFac%s_%sTesla_relaxationRatesScaledUp', ...
    whichLipid,r1MatlabSimulationDate,scalFacMatlabSimulationDate, ...
    fieldstrengthString);
scaledUpR1_path = [savingPath savingName '.mat'];
save(scaledUpR1_path,'r1ScaledUp_theta_phi_NN', ...
    'r1ScaledUp_theta_NN','r1ScaledUp_NN','whichLipid', ...
    'fieldStrength','scalFacMatlabSimulationDate', ...
    'r1MatlabSimulationDate','nearestNeighbourCases', ...
    'orientationAngles','positionAngles','r1GromascSimulationDate', ...
    'scalFacGromacsSimulationDate');
end
