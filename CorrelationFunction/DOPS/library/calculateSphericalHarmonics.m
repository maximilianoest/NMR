function [F1,F2] = calculateSphericalHarmonics(theta,phi ...
    ,nearestNeighbourDistancesPow3)
% This function calculates the spherical harmonics based on the angles and
% distances.

SIN = sin(theta);
COS = cos(theta);

% spherical harmonics
F1 = SIN.*COS.*exp(-1i*phi)./nearestNeighbourDistancesPow3;
F2 = SIN.^2.*exp(-2i*phi)./nearestNeighbourDistancesPow3;
end
