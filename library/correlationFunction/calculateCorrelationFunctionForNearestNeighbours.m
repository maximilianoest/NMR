function [correlationFunction] = ...
    calculateCorrelationFunctionForNearestNeighbours( ...
    sphericalHarmonic,numLags)
% EXPLANATIONS OF CALCULATION:
% IFFt is slower than FFT, because FFT has a build in function which
% detects whether the input data is real or complex and therefore chooses
% an optimized way to calculate the inverse FT. Because the input data of
% IFFT is real, this data can be seen as conjugate symmetric. Therefore,
% the 'symmetric' flag at IFFT is set.
% The reason why in frequency domain the squared real and imaginary part
% are summed is, because the autocorrelation require the "convolution" (not
% directly the convolution, but more the correlation) of a signal and its 
% complex conjugate. This is a multiplication of the signal
% and its complex conjugate in frequency domain and therefore the real and
% imaginary part summed up.


[~,timeSteps] = size(sphericalHarmonic);
zeroPaddingLength = 2^(nextpow2(timeSteps));

fftSphericalHarmonic = fft(sphericalHarmonic,zeroPaddingLength,2);

multiplicationFunction = real(fftSphericalHarmonic).^2 ...
    + imag(fftSphericalHarmonic).^2;
clearvars fftSphericalHarmonic

correlationFunction = ifft(multiplicationFunction,[],2,'symmetric');
correlationFunction = correlationFunction(:,1:numLags)/timeSteps;

end

