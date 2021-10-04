function [sumCorrelationFunction] = ...
    calculateCorrelationFunction(sphericalHarmonic,numLags)

tic
[~,timeSteps] = size(sphericalHarmonic);
zeroPaddingLength = 2^(nextpow2(timeSteps)+1);

fftSphericalHarmonic = fft(cast(sphericalHarmonic,'single') ...
    ,zeroPaddingLength,2);

multiplicationFunction = real(fftSphericalHarmonic).^2 ...
    + imag(fftSphericalHarmonic).^2;

correlationFunction = ifft(multiplicationFunction,[],2);

sumCorrelationFunction = sum(correlationFunction(:,1:numLags),1)/timeSteps;

% fft creates a vector with higher amplitudes => to get a valid result from
% fft you have to devide by the length of the vector
% ifft reduces to the original amplitudes i.e. fft(x) = higher amplitudes ,
% ifft(fft(x)) = original amplitudes => with ifft(fft(x) *fft(y)) you need
% to devide by the length but only once
toc
 end



