function scalFacSaving_path = determineScalingFactors( ...
    r1ForDifferentNN_path,fieldStrength)

%% initialize system
% Set Up parameters

results = load(r1ForDifferentNN_path);

% calculate R1 from correlation function. Here will be implemented a
% function which allows to calculate R1 directly from correlation function.
% For now, the R1s will be taken directly from the results files.

% simulation information
whichLipid = results.whichLipid;
matlabSimulationDate = results.matlabSimulationDate;
gromacsSimulationDate = results.gromacsSimulationDate;
orientationAngles = results.orientationAngles;
positionAngles = results.positionAngles;
nearestNeighbourCases = results.nearestNeighbourCases;
atomCounter = results.atomCounter;
gromacsFileName = results.gromacsFileName;

% unscaled r1s calculated from the correlation functions.
r1_theta_phi_NN = results.r1_theta_phi_NN;
r1_theta_NN = squeeze(mean(r1_theta_phi_NN,2));
r1_NN = squeeze(mean(r1_theta_NN,1));

%% determine scaling factors

scalFac_theta_phi_NN = r1_theta_phi_NN(:,:,1)./r1_theta_phi_NN;
scalFac_theta_NN = r1_theta_NN(:,1)./r1_theta_NN;
scalFac_NN = r1_NN(1)./r1_NN;

% validation to make sure the scaling factors are right.
upscaledR1_theta_phi_NN = scalFac_theta_phi_NN.*r1_theta_phi_NN;
upscaledR1_theta_NN = scalFac_theta_NN.*r1_theta_NN;
upscaledR1_NN = scalFac_NN.*r1_NN;

checkSum_theta_phi_NN = sum(sum(sum(abs(upscaledR1_theta_phi_NN  ...
    - squeeze(r1_theta_phi_NN(:,:,1)))))); 
checkSum_theta_NN = sum(sum(abs(upscaledR1_theta_NN ...
    - squeeze(r1_theta_NN(:,1)))));
checkSum_NN = sum(abs(upscaledR1_NN-upscaledR1_NN(1)));

if checkSum_theta_phi_NN > 1e-10 || checkSum_theta_NN > 1e-10 ...
        || checkSum_NN > 1e-10
    warning('determineScalingFactors:checkSumTooLarge', ...
        'Warning! Something went wrong when scaling factors were calculated');
else
    disp('All check sums for scaling factors are fine');
    savingPath = initializeSystemForSavingScalingFactors();
    fieldstrengthString = strrep(num2str(fieldStrength),'.','');
    savingName = sprintf('%s_r1%s_%sTesla_scalingFactors', ...
        whichLipid,matlabSimulationDate,fieldstrengthString);
    scalFacSaving_path = [savingPath savingName];
    save(scalFacSaving_path,'scalFac_theta_phi_NN', ...
        'scalFac_theta_NN','scalFac_NN','fieldStrength', ...
        'whichLipid','matlabSimulationDate', ...
        'orientationAngles','positionAngles','nearestNeighbourCases', ...
        'gromacsSimulationDate','atomCounter','gromacsFileName');
   
end

end







