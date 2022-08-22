function testScalingFactorsAndR1sToScaleForCompatibility( ...
    scalFacResults,r1ForDifferentNN)

% check which lipids are used
scalFacWichLipid = scalFacResults.whichLipid;
r1WhichLipid = r1ForDifferentNN.whichLipid;
if ~strcmp(scalFacWichLipid,r1WhichLipid)
    error('scaleUpR1WithSCalingFactors:differentLipid', ...
        ['The lipid which determined the scaling factors is not the same' ...
        ' like the one for which R1 should be upscaled']);
end

% check which field strengths are used
scalFacFieldStrength = scalFacResults.fieldStrength;
r1FieldStrength = r1ForDifferentNN.fieldStrength;
if ~(scalFacFieldStrength == r1FieldStrength)
    error('scaleUpR1WithSCalingFactors:differntFieldstrength', ...
        ['The fieldstrength from determining the scaling factor differs' ...
        'from the fieldstrength for which the R1s are calculated']);
end

% check which simulations are used
scalFacStartDateOfSimulation = scalFacResults.matlabSimulationDate;
r1StartDateOfSimulation = r1ForDifferentNN.matlabSimulationDate; 
if strcmp(scalFacStartDateOfSimulation,r1StartDateOfSimulation)
    warning('scaleUpR1WithSCalingFactors:sameDateWarning', ...
        ['The scaling factors are from a simulation with the same date' ...
        ' like the r1 rates that should be upscaled']);
end

% % check for different NN cases
% scalFacNearestNeighbourCases = scalFacResults.nearestNeighbourCases;
% r1NearestNeighbourCases = r1ForDifferentNN.nearestNeighbourCases;

% % check whether the orientation angles are the same
% scalFacOrientationAngles = scalFacResults.orientationAngles;
% r1OrientationAngles = r1ForDifferentNN.orientationAngles;
% for oriNr = 1:length(scalFacOrientationAngles)
%     if scalFacOrientationAngles(oriNr) ~= r1OrientationAngles(oriNR)
%         error('testingCompatibilityOfScalFac:unequalOrientations', ...
%             ['The orientation angle for 
%     end
% end
% 
% % check whether the position angles are the same
% 
% 
% 
% % check whehter the nearest neighbour cases are the same
% 
% 
% scalFacPositionAngles = scalFacResults.positionAngles;
% 
% 
% 
% 
% 
% 
% r1PositionAngles = r1ForDifferentNN.positionAngles;
% 
% 
% 

end
