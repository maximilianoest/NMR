function [F1,F2] = calculateSphericalHarmonics(nearestNeighboursX ...
    ,nearestNeighboursY,nearestNeighboursZ,nearestNeighbourDistancesPow3)
% This function first calculates the transformation from cartesian to 
% spherical coordinates. Afterwards the spherical harmonics are calculated.

% tranform cartesian to spherical coordinates
hypotuseXY = hypot(nearestNeighboursX,nearestNeighboursY);
theta = pi/2-atan2(nearestNeighboursZ,hypotuseXY);
phi = atan2(nearestNeighboursY,nearestNeighboursX);

SIN = sin(theta);
COS = cos(theta);

% spherical harmonics
F1 = SIN.*COS.*exp(-1i*phi)./nearestNeighbourDistancesPow3;
F2 = SIN.^2.*exp(-2i*phi)./nearestNeighbourDistancesPow3;
end
