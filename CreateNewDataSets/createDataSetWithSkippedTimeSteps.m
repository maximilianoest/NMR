clc
clear all

%% Lipid Data
disp('Started loading lipid data')
data = load(['/home/fschyboll/Dissertation/Paper4/Gromacs_Final' ...
    '/prd/Traj_Extracted/Lipid_H_500ns_1ps_wh.mat']);
path2Save = '/home/moesterreicher/Data/Myelin/';
disp('Finished loading lipid data')

dataSet = data.Lipid_H_500ns_1ps_wh;
clearvars -except dataSet path2Save
sizeOfData = size(dataSet);
disp(sizeOfData)

newStepSize = 2;
Lipid_H_500ns_2ps_wh = dataSet(:,:,1:newStepSize:end);
dataContentName = 'Lipid_H_500ns_2ps_wh';
save([path2Save dataContentName '.mat'],'Lipid_H_500ns_2ps_wh','-v7.3');
clearvars -except dataSet path2Save

newStepSize = 3;
Lipid_H_500ns_3ps_wh = dataSet(:,:,1:newStepSize:end);
dataContentName = 'Lipid_H_500ns_3ps_wh';
save([path2Save dataContentName '.mat'],'Lipid_H_500ns_3ps_wh','-v7.3');
clearvars -except dataSet path2Save

newStepSize = 4;
Lipid_H_500ns_4ps_wh = dataSet(:,:,1:newStepSize:end);
dataContentName = 'Lipid_H_500ns_4ps_wh';
save([path2Save dataContentName '.mat'],'Lipid_H_500ns_4ps_wh','-v7.3');
clearvars

%% Water data
disp('Started loading water data')
data = load(['/home/fschyboll/Dissertation/Paper3/Matlab' ...
        '/Eval_MonoLayer3/Traj_mat/Short_50ns_025ps' ...
        '/water_H_50ns_025ps_wh.mat']);
path2Save = '/home/moesterreicher/Data/Myelin/';
disp('Finished loading water data')

dataSet = data.water_H_50ns_025ps_wh;
clearvars -except dataSet path2Save
sizeOfData = size(dataSet);
disp(sizeOfData)

newStepSize = 2;
water_H_50ns_05ps_wh = dataSet(:,:,1:newStepSize:end);
dataContentName = 'water_H_50ns_05ps_wh';
save([path2Save dataContentName '.mat'],'water_H_50ns_05ps_wh','-v7.3');
clearvars -except dataSet path2Save

newStepSize = 3;
water_H_50ns_075ps_wh = dataSet(:,:,1:newStepSize:end);
dataContentName = 'water_H_50ns_075ps_wh';
save([path2Save dataContentName '.mat'],'water_H_50ns_075ps_wh','-v7.3');
clearvars -except dataSet path2Save

newStepSize = 4;
water_H_50ns_1ps_wh = dataSet(:,:,1:newStepSize:end);
dataContentName = 'water_H_50ns_1ps_wh';
save([path2Save dataContentName '.mat'],'water_H_50ns_1ps_wh','-v7.3');
clearvars -except dataSet path2Save
