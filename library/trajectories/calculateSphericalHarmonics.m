function [zerothOrderSphericalHarmonic,firstOrderSphericalHarmonic ...
    ,secondOrderSphericalHarmonic] = calculateSphericalHarmonics( ...
    polarAngle,azimuthAngle,nearestNeighbourDistancesPow3)

zerothOrderSphericalHarmonic = 3/4*(1-3*(cos(polarAngle)).^2)./ ...
    nearestNeighbourDistancesPow3;
firstOrderSphericalHarmonic = sin(polarAngle).*cos(polarAngle) ...
    .*exp(-1i*azimuthAngle)./nearestNeighbourDistancesPow3;
secondOrderSphericalHarmonic = sin(polarAngle).^2 ...
    .*exp(-2i*azimuthAngle)./nearestNeighbourDistancesPow3;
end
