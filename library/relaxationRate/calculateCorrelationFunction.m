function [correlationFunction] = ...
    calculateCorrelationFunction(sphericalHarmonic,numLags)

[~,timeSteps] = size(sphericalHarmonic);
zeroPaddingLength = 2^(nextpow2(timeSteps)+1);

% Convolution
tic 
correlationFunction = ifft(fft(sphericalHarmonic,zeroPaddingLength,2) ... 
    .*conj(fft(sphericalHarmonic,zeroPaddingLength,2)),[],2)/timeSteps;
correlationFunction = sum(correlationFunction(:,1:numLags),1);
toc 


% fft creates a vector with higher amplitudes => to get a valid result from
% fft you have to devide by the length of the vector
% ifft reduces to the original amplitudes i.e. fft(x) = higher amplitudes ,
% ifft(fft(x)) = original amplitudes => with ifft(fft(x) *fft(y)) you need
% to devide by the length but only once

end



