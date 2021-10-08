clc
clear all

load('C:\Users\maxoe\Google Drive\Promotion\Results\PerformanceAnalysis\20211001_Results_PerformanceAnalysis_Lipid_H_500ns_1ps_nH40.mat')
oldR1 = r1WithPerturbationTheory;
oldTimeTracks = timeTracks;
oldAtomCounter = atomCounter;


load('C:\Users\maxoe\Google Drive\Promotion\Results\PerformanceAnalysis\20211008_Results_PerformanceAnalysis_Lipid_H_500ns_1ps_nH40.mat')
newR1 = r1WithPerturbationTheory;
newTimeTracks = timeTracks;
newAtomCounter = atomCounter;

