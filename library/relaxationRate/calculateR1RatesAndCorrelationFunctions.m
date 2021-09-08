function [r1WithPerturbationTheory,correlationFunction1W0Saver ...
    ,correlationFunction2W0Saver] = ...
    calculateR1RatesAndCorrelationFunctions(rotatedX,rotatedY,rotatedZ ...
    ,nearestNeighbourDistancesPow3,lags,omega0,deltaT,DD ...
    ,shiftForCorrelationFunction)

polarAngle = zeros(size(rotatedX,1),size(rotatedY,2));
azimuthAngle = zeros(size(rotatedX,1),size(rotatedY,2));

firstOrderSphericalHarmonic = zeros(size(rotatedX,1),size(rotatedY,2));
secondOrderSphericalHarmonic = zeros(size(rotatedX,1),size(rotatedY,2));
[nearestNeighbours,timeSteps] = size(firstOrderSphericalHarmonic);
zeroPaddingLength = 2^(nextpow2(timeSteps)+1);

correlationFunction1W0 = zeros(nearestNeighbours,zeroPaddingLength);
correlationFunction2W0 = zeros(nearestNeighbours,zeroPaddingLength);

spectralDensity1W0 = single(0);
spectralDensity2W0 = single(0);
r1WithPerturbationTheory = single(0);


transformToSphericalCoordinates();
    function transformToSphericalCoordinates()
        
        hypotuseXY = hypot(rotatedX,rotatedY);
        polarAngle = pi/2-atan2(rotatedZ,hypotuseXY);
        azimuthAngle = atan2(rotatedY,rotatedX);
        
    end

calculateSphericalHarmonics();
    function calculateSphericalHarmonics()
        
        firstOrderSphericalHarmonic = sin(polarAngle).*cos(polarAngle) ...
            .*exp(-1i*azimuthAngle)./nearestNeighbourDistancesPow3;
        secondOrderSphericalHarmonic = sin(polarAngle).^2 ...
            .*exp(-2i*azimuthAngle)./nearestNeighbourDistancesPow3;
    end

calculateCrossCorrelationFunction();
    function calculateCrossCorrelationFunction()
        
        firstOrderSphericalHarmonic = fft(firstOrderSphericalHarmonic ...
            ,zeroPaddingLength,2);
        % Convolution
        correlationFunction1W0 = ifft(firstOrderSphericalHarmonic ...
            .*conj(firstOrderSphericalHarmonic),[],2)/timeSteps;
        correlationFunction1W0 = sum(correlationFunction1W0(:,1:lags));
        
        secondOrderSphericalHarmonic = fft(secondOrderSphericalHarmonic ...
            ,zeroPaddingLength,2);
        % Convolution
        correlationFunction2W0 = ifft(secondOrderSphericalHarmonic ...
            .*conj(secondOrderSphericalHarmonic),[],2)/timeSteps;
        correlationFunction2W0 = sum(correlationFunction2W0(:,1:lags));
    end

correlationFunction1W0Saver = single(correlationFunction1W0( ...
    1:shiftForCorrelationFunction:end));
correlationFunction2W0Saver = single(correlationFunction2W0( ...
    1:shiftForCorrelationFunction:end));

calculateSpectralDensities();
    function calculateSpectralDensities()
        
        spectralDensity1W0 = 2*(deltaT*cumsum(correlationFunction1W0 ...
            .*exp(-1i*omega0*deltaT*(0:lags-1))));
        spectralDensity2W0 = 2*(deltaT*cumsum(correlationFunction2W0 ...
            .*exp(-2i*omega0*deltaT*(0:lags-1))));
        
        spectralDensity1W0 = spectralDensity1W0(end);
        spectralDensity2W0 = spectralDensity2W0(end);
        
    end

calculateR1WithSpectralDensity();
    function calculateR1WithSpectralDensity()
        
        r1WithPerturbationTheory = DD*3/2*( ...
            abs(real(mean(spectralDensity1W0))) ...
            +abs(real(mean(spectralDensity2W0))));
        
    end


end
