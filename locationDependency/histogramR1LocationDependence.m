clc; clear all; close all;

randomRadii = rand(1,10000000)+1;
randomR1 = 1./randomRadii.^6;

histogram(randomR1,'BinWidth',0.01);

