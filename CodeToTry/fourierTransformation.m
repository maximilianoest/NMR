clc
clear all

samplingFrequenz = 100;
dt = 1/samplingFrequenz;
timeAxis = 0:dt:5;
frequenz = 30;
signal = sin(2*pi*frequenz*timeAxis);

frequenzAchse = 0:
fftSignal = fftshift(fft(signal));

figure(1)
plot(timeAxis,signal)
grid on

figure(2)
plot(real(fftSignal(end/2:end)))
grid on