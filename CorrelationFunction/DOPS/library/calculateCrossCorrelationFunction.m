function [crosscorrelationFunction] = ...
    calculateCrossCorrelationFunction(sphericalHarmonic ... 
    ,complexConjugatedSphericalHarmonic,numLags)
% This function calculates the cross correlation function of the spherical
% harmonics.

[~,timeSteps] = size(sphericalHarmonic);
zeroPaddingLength = 2^(nextpow2(timeSteps)+1);

firstSphericalHarmonicFFT = fft(sphericalHarmonic,zeroPaddingLength,2); 
secondSphericalHarmonicFFT = fft( ...
    conj(complexConjugatedSphericalHarmonic),zeroPaddingLength,2); 

crosscorrelationFunction = ifft(firstSphericalHarmonicFFT ... 
    .*conj(secondSphericalHarmonicFFT),[],2);
crosscorrelationFunction = crosscorrelationFunction(:,(1:numLags+1) ) ...
    /timeSteps;
end



