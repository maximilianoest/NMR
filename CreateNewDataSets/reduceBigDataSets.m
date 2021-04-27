clc
clear all %#ok<CLALL>

%% load lipid Data
dataLoaded = false;
if ~dataLoaded
    data = load(['/home/fschyboll/Dissertation/Paper4/Gromacs_Final' ...
        '/prd/Traj_Extracted/Lipid_H_500ns_1ps_wh.mat']);
end

%% analyse data
dataSet = data.Lipid_H_500ns_1ps_wh;
clearvars -except dataSet
sizeOfData = size(dataSet);
disp(sizeOfData)
numberOfAtoms = sizeOfData(1);

%% extract small data set to test on local machine
numberOfAtomsInSmallDataSet = 40;
randomStartingPointSmallDataSet = floor((numberOfAtoms ...
    -numberOfAtomsInSmallDataSet)*rand());
Lipid_H_500ns_1ps_nH40 = dataSet( ...
    randomStartingPointSmallDataSet:randomStartingPointSmallDataSet ...
    +numberOfAtomsInSmallDataSet,:,:);

%% extract medium data set for calculations on local machine/server
numberOfAtomsInMediumDataSet = 200;
randomStartingPointMediumDataSet = floor((numberOfAtoms ...
    -numberOfAtomsInMediumDataSet)*rand());
Lipid_H_500ns_1ps_nH200 = dataSet( ...
    randomStartingPointMediumDataSet:randomStartingPointMediumDataSet ...
    +numberOfAtomsInMediumDataSet,:,:);

%% extract big data set for calculations on server
numberOfAtomsInBigDataSet = 1000;
randomStartingPointBigDataSet = floor((numberOfAtoms ...
    -numberOfAtomsInBigDataSet)*rand());
Lipid_H_500ns_1ps_nH1000 = dataSet( ...
    randomStartingPointBigDataSet:randomStartingPointBigDataSet ...
    +numberOfAtomsInBigDataSet,:,:);

%% save data sets
path2Save = '/home/moesterreicher/Data/Myelin/';
dataContentName = 'Lipid_H_500ns_1ps_nH';
save([path2Save dataContentName num2str( ...
    numberOfAtomsInSmallDataSet) '.mat'],'Lipid_H_500ns_1ps_nH40' ...
    ,'randomStartingPointSmallDataSet','-v7.3');

save([path2Save dataContentName num2str( ...
    numberOfAtomsInMediumDataSet) '.mat'],'Lipid_H_500ns_1ps_nH200' ...
    ,'randomStartingPointMediumDataSet','-v7.3');

save([path2Save dataContentName num2str( ...
    numberOfAtomsInBigDataSet) '.mat'],'Lipid_H_500ns_1ps_nH1000' ...
    ,'randomStartingPointBigDataSet','-v7.3');

clc
clear all %#ok<CLALL>

%% load water data
dataLoaded = false;
if ~dataLoaded
    data = load(['/home/fschyboll/Dissertation/Paper3/Matlab' ...
        '/Eval_MonoLayer3/Traj_mat/Short_50ns_025ps' ...
        '/water_H_50ns_025ps_wh.mat']);
end

%% analyse data
dataSet = data.water_H_50ns_025ps_wh;
clearvars -except dataSet
sizeOfData = size(dataSet);
disp(sizeOfData)
numberOfAtoms = sizeOfData(1);

%% extract small data set to test on local machine
numberOfAtomsInSmallDataSet = 100;
randomStartingPointSmallDataSet = floor( ...
    (numberOfAtoms-numberOfAtomsInSmallDataSet)*rand());
if mod(randomStartingPointSmallDataSet,2) ~= 0
   randomStartingPointSmallDataSet = randomStartingPointSmallDataSet-1; 
end
Water_H_50ns_025ps_nH100 = dataSet(randomStartingPointSmallDataSet ...
    :randomStartingPointSmallDataSet+numberOfAtomsInSmallDataSet,:,:);

%% extract medium data set for calculations on local machine/server
numberOfAtomsInMediumDataSet = 500;
randomStartingPointMediumDataSet = floor( ...
    (numberOfAtoms-numberOfAtomsInMediumDataSet)*rand());
if mod(randomStartingPointMediumDataSet,2) ~= 0
   randomStartingPointMediumDataSet = randomStartingPointMediumDataSet-1;
end
Water_H_50ns_025ps_nH500 = dataSet(randomStartingPointMediumDataSet ...
    :randomStartingPointMediumDataSet+numberOfAtomsInMediumDataSet,:,:);

%% extract big data set for calculations on server
numberOfAtomsInBigDataSet = 2500;
randomStartingPointBigDataSet = floor( ...
    (numberOfAtoms-numberOfAtomsInBigDataSet)*rand());
if mod(randomStartingPointBigDataSet,2) ~= 0
   randomStartingPointBigDataSet = randomStartingPointBigDataSet-1; 
end
Water_H_50ns_025ps_nH2500 = dataSet(randomStartingPointBigDataSet ...
    :randomStartingPointBigDataSet+numberOfAtomsInBigDataSet,:,:);

%% save data sets
path2Save = '/home/moesterreicher/Data/Myelin/';
dataContentName = 'Water_H_50ns_025ps_nH';
save([path2Save dataContentName num2str( ...
    numberOfAtomsInSmallDataSet) '.mat'],'Water_H_50ns_025ps_nH100' ...
    ,'randomStartingPointSmallDataSet','-v7.3');
save([path2Save dataContentName num2str( ...
    numberOfAtomsInMediumDataSet) '.mat'],'Water_H_50ns_025ps_nH500' ...
    ,'randomStartingPointMediumDataSet','-v7.3');
save([path2Save dataContentName num2str( ...
    numberOfAtomsInBigDataSet) '.mat'],'Water_H_50ns_025ps_nH2500' ...
    ,'randomStartingPointBigDataSet','-v7.3');





