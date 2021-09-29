function [correlationFunction] = ...
    calculateCorrelationFunction(sphericalHarmonic,numLags)

[~,timeSteps] = size(sphericalHarmonic);
zeroPaddingLength = 2^(nextpow2(timeSteps)+1);

% Convolution
tic 
correlationFunction = ifft(fft(sphericalHarmonic,zeroPaddingLength,2) ... 
    .*fft(conj(sphericalHarmonic),zeroPaddingLength,2),[],2)/timeSteps;
correlationFunction = sum(correlationFunction(:,1:numLags),1);
toc 


end



