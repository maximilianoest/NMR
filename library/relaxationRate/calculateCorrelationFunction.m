function [correlationFunction] = ...
    calculateCorrelationFunction(sphericalHarmonic,numLags)

[~,timeSteps] = size(sphericalHarmonic);
zeroPaddingLength = 2^(nextpow2(timeSteps)+1);

% Convolution
correlationFunction = (ifft(fft(sphericalHarmonic,zeroPaddingLength,2) ... 
    .*fft(conj(sphericalHarmonic),zeroPaddingLength,2),[],2)/timeSteps);
correlationFunction = correlationFunction(:,1:numLags);

end



