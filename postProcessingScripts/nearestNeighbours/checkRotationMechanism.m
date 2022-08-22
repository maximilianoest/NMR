clc;clear all;close all;%#ok<*NOPTS>
configuration = readConfigurationFile('nNConfig.txt');
results = load(configuration.corrFuncDOPSForScalFac_path);

orAng = rad2deg(results.orientationAngle)
posAng = rad2deg(results.positionAngle) 

invRotOr = results.inverseRotationMatrixOrientation
invRotPos = results.inverseRotationMatrixPosition

rotOr = results.rotationMatrixOrientation
rotPos = results.rotationMatrixPosition

totalRot = rotOr*rotPos
vecX = [1 0 0];
vecY = [0 1 0];
vecZ = [0 0 1];

rotVecX = totalRot*vecX'
rotVecY = totalRot*vecY'
rotVecZ = totalRot*vecZ'

zeroRot = rotOr*rotPos*invRotPos*invRotOr
