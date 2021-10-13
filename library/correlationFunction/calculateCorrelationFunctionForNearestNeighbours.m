function [correlationFunction] = ...
    calculateCorrelationFunctionForNearestNeighbours( ...
    sphericalHarmonic,numLags)

[~,timeSteps] = size(sphericalHarmonic);
zeroPaddingLength = 2^(nextpow2(timeSteps)+1);

fftSphericalHarmonic = fft(sphericalHarmonic,zeroPaddingLength,2);

multiplicationFunction = real(fftSphericalHarmonic).^2 ...
    + imag(fftSphericalHarmonic).^2;

correlationFunction = ifft(multiplicationFunction,[],2);
correlationFunction = correlationFunction(:,1:numLags)/timeSteps;

end

