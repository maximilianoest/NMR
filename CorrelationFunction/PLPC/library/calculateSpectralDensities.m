function [w1,w2] = calculateSpectralDensities(correlationFunction1, ...
    correlationFunction2,omega0,deltaT,lags)
% This function calculates the spectral densities of the given
% correlation functions

w1 = 2*(deltaT*cumsum(sum(correlationFunction1) ... 
        .*exp(-1i*omega0*deltaT*(0:lags-1))));
w2 = 2*(deltaT*cumsum(sum(correlationFunction2) ...
        .*exp(-2i*omega0*deltaT*(0:lags-1))))  ;

end
