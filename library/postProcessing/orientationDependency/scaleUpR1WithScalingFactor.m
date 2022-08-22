function scaledUpR1_path = scaleUpR1WithScalingFactor(unscaledR1_path, ...
    scalFac_path)
r1Data = load(unscaledR1_path);
scalFacData = load(scalFac_path);
fprintf('%25s:       R1   | scaling Factors \n','SIMULATION INFORMATION');
fprintf('---------------------------------------------------------\n');

r1WhichLipid = r1Data.whichLipid;
scalFacWhichLipid = scalFacData.whichLipid;
fprintf('  %23s: %8s   | %8s\n','Lipid',r1WhichLipid,scalFacWhichLipid);

r1NN = r1Data.nearestNeighbours;
scalFacNNCases = scalFacData.nearestNeighbourCases;
scalFacIndex = find(r1NN == scalFacNNCases);
fprintf('  %23s: %8i   | %8i (index: %i)\n','Nearest neighbours',r1NN, ...
    scalFacNNCases(scalFacIndex),scalFacIndex);

r1MatlabSimDate = r1Data.matlabSimulationDate;
scalFacMatlabSimDate = scalFacData.matlabSimulationDate;
fprintf('  %23s: %s   | %s\n','Matlab simulation date',r1MatlabSimDate, ...
    scalFacMatlabSimDate);

r1GromascSimDate = r1Data.gromacsSimulationDate;
scalFacGromacsSimDate = r1Data.gromacsSimulationDate;
fprintf('  %23s: %s   | %s\n','Gromacs simulation date', ...
    r1GromascSimDate,scalFacGromacsSimDate);

r1FieldStrength = r1Data.fieldStrength;
scalFacFieldStrength = scalFacData.fieldStrength;
fprintf('  %23s: %8.2f   | %8.2f\n','Field strength [T]',r1FieldStrength, ...
    scalFacFieldStrength);

r1AtomCounter = r1Data.atomCounter;
scalFacAtomCounter = scalFacData.atomCounter;
fprintf('  %23s: %8i   | %8i\n','Calculated atoms',r1AtomCounter, ...
    scalFacAtomCounter);

fprintf('---------------------------------------------------------\n');
scalingFactor = scalFacData.scalFac_NN(scalFacIndex);
fprintf(' SCALING FACTOR: %f\n',scalFacData.scalFac_NN(scalFacIndex));

scaledUpR1_theta_phi = r1Data.unscaledR1_theta_phi*scalingFactor;
r1Data.scaledUpR1_theta_phi = scaledUpR1_theta_phi;
r1Data.scalingFactorInformation = scalFacData;
r1Data.usedScalingFactor = scalingFactor;

savingPath = initializeSystemForSavingR1();
fieldstrengthString = strrep(num2str(r1FieldStrength),'.','');
savingName = sprintf('%s_%sTesla_%iH_%iNN_r1%s_scaledRelaxationRatesR1', ...
    r1WhichLipid,fieldstrengthString,r1AtomCounter,r1NN, ...
    r1MatlabSimDate);
scaledUpR1_path = [savingPath savingName '.mat'];
save(scaledUpR1_path,'-struct','r1Data');

end
