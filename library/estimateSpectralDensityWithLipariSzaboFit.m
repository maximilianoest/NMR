function [spectralDensityLiSz] = ...
    estimateSpectralDensityWithLipariSzaboFit(correlationFunction,omega ...
    ,deltaT,outputLogFileName)

[numberOfHs,lags] = size(correlationFunction);
timeAxis = (0:lags-1)';

spectralDensitysLiSz = zeros(1,numberOfHs);

firstEstimation = @(b,x) b(1)^2*exp(-x/b(3))+b(2);
secondEstimation = @(b,x) (1-b(1)^2)*exp(-x/b(2))+b(1)^2*exp(-x/b(3));
fileId = fopen(outputLogFileName,'a');

opts = optimset('Display','off');
for atomNumber = 1:numberOfHs
    firstFitLimit = 0.004; %1000*picoSecond/(timeSteps*deltaT); %3500
    secondFitLimit = 0.8;
    normalizedCurve = (correlationFunction(atomNumber,:)') ...
        /correlationFunction(atomNumber,1);
    
    while true
        try
            firstEstimationParameters = lsqcurvefit(firstEstimation ...
                ,[1,200,0],timeAxis(1:1:round(firstFitLimit*lags)) ...
                ,normalizedCurve(1:1:round(firstFitLimit*lags)) ...
                ,[0 0 20],[1 1 600],opts);
            break
        catch
            fprintf(fileId ...
                ,['Error while fitting first estimation curve. ' ...
                ,'Trying with lower first limit.\n']);
            firstFitLimit = firstFitLimit-0.0002;
        end
    end
    
    while true
        try
            secondEstimationParameters = lsqcurvefit(secondEstimation ...
                ,[0.5,5000,firstEstimationParameters(3)] ...
                ,timeAxis(1:4:round(secondFitLimit*lags)) ...
                ,normalizedCurve(1:4:round(secondFitLimit*lags)) ...
                ,[0 100 0.5*firstEstimationParameters(3)] ...
                ,[1 2000000 1.5*firstEstimationParameters(3)],opts);
            break
        catch
            fprintf(fileId ...
                ,['Error while fitting first estimation curve. ' ...
                ,'Trying with lower second limit.\n']);
            secondFitLimit = secondFitLimit-0.05;
        end
    end
    
    orderParameterSQRT = secondEstimationParameters(1);
    internalCorrelationTime = secondEstimationParameters(2)*deltaT;
    globalCorrelationTime = secondEstimationParameters(3)*deltaT;
    
    spectralDensitysLiSz(atomNumber) = 2*correlationFunction(atomNumber,1) ...
        *(orderParameterSQRT^2*globalCorrelationTime ...
        /(1+(globalCorrelationTime*omega)^2) ...
        +(1-orderParameterSQRT^2)*internalCorrelationTime ...
        /(1+(internalCorrelationTime*omega)^2));
    
%     firstFit = firstEstimationParameters(1)^2 ...
%         *exp(-timeAxis(1:round(0.1*lags))/firstEstimationParameters(3)) ...
%         +firstEstimationParameters(2);
%     secondFit = (1-secondEstimationParameters(1)^2) ...
%         *exp(-timeAxis/secondEstimationParameters(2)) ...
%         +secondEstimationParameters(1)^2 ...
%         *exp(-timeAxis/secondEstimationParameters(3));
%     
%     plot(real(secondFit))
%     hold on
%     plot(real(normalizedCurve))
%     hold off
%     pause(0.2)
%     
%     plot(abs(real(firstFit)))
%     hold on
%     plot(abs(real(normalizedCurve(1:round(0.1*lags)))),'x')
%     hold off
%     pause(0.2)

end

spectralDensityLiSz = sum(spectralDensitysLiSz);
fclose(fileId);

end
