function [correlationFunction] = ...
    calculateCorrelationFunctionForNearestNeighbours( ...
    sphericalHarmonic,numLags)

[~,timeSteps] = size(sphericalHarmonic);
zeroPaddingLength = 2^(nextpow2(timeSteps)+1);

correlationFunction = ifft(fft(sphericalHarmonic,zeroPaddingLength,2) ... 
    .*conj(fft(sphericalHarmonic,zeroPaddingLength,2)),[],2)/timeSteps;
correlationFunction = sum(correlationFunction(:,1:numLags),1);

end

