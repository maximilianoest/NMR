function [correlationFunction] = ...
    calculateCorrelationFunction(sphericalHarmonic,numLags)

[~,timeSteps] = size(sphericalHarmonic);
zeroPaddingLength = 2^(nextpow2(timeSteps)+1);

% Convolution
logMemoryUsage('C:\Users\maxoe\Google Drive\Promotion\Results\OrientationDependency\20210927_LogFile_OrientationDependency_Lipid_H_500ns_1ps_nH40.txt')
correlationFunction = (ifft(fft(sphericalHarmonic,zeroPaddingLength,2) ... 
    .*fft(conj(sphericalHarmonic),zeroPaddingLength,2),[],2)/timeSteps);
logMemoryUsage('C:\Users\maxoe\Google Drive\Promotion\Results\OrientationDependency\20210927_LogFile_OrientationDependency_Lipid_H_500ns_1ps_nH40.txt')
correlationFunction = sum(correlationFunction(:,(1:numLags)),1);

logMemoryUsage('C:\Users\maxoe\Google Drive\Promotion\Results\OrientationDependency\20210927_LogFile_OrientationDependency_Lipid_H_500ns_1ps_nH40.txt')
end



