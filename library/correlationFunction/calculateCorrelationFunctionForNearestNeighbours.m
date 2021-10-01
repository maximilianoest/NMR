function [correlationFunction] = ...
    calculateCorrelationFunctionForNearestNeighbours( ...
    sphericalHarmonic,numLags)

[~,timeSteps] = size(sphericalHarmonic);
zeroPaddingLength = 2^(nextpow2(timeSteps)+1);

fftSphericalHarmonic = fft(sphericalHarmonic,zeroPaddingLength,2);

correlationFunction = ifft(fftSphericalHarmonic.*conj( ...
    fftSphericalHarmonic),[],2);
correlationFunction = correlationFunction(:,1:numLags)/timeSteps;

end

