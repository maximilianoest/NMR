clc
clear all %#ok<CLALL>
close all

data = load(['C:\Users\maxoe\Google Drive\Promotion\Data' ...
    '\Myelin\Lipid_RelaxationRates_nH40.mat']);


% TODO: wirte function getFieldNamesOfStruct(structName) -> put it in 
% library in subfolder (see if it works)
fieldNames = fieldnames(data);
if length(fieldNames) == 1
    data = data.(fieldNames{1});
    fieldNames = fieldnames(data);
end

r1Rate = data.r1WithPerturbationTheory;




 



