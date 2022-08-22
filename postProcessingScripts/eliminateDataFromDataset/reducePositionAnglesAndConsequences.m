clc;close all;clear all;
PLPC_path = "C:\Users\maxoe\Google Drive\Promotion\Results\LIPIDS\PLPC\orientationDependency\20220502_Results_orientationDependency_PLPClipid.mat";
PSM_path = "C:\Users\maxoe\Google Drive\Promotion\Results\LIPIDS\PSM\orientationDependency\20220530_Results_orientationDependency_PSMlipid.mat";
DOPS_path = "C:\Users\maxoe\Google Drive\Promotion\Results\LIPIDS\DOPS\orientationDependency\20220713_Results_orientationDependency_DOPSlipid.mat";

createCoincidingResultFiles(DOPS_path,PLPC_path);
createCoincidingResultFiles(DOPS_path,PSM_path);

function createCoincidingResultFiles(reference_path,compare_path)
referenceData = load(reference_path);
referencePhiAngles = referenceData.positionAngles;

compareData = load(compare_path);
comparePhiAngles = compareData.positionAngles;

coincidentIndex = [];


for refIndex = 1:numel(referencePhiAngles) 
    
    for compareIndex = 1:numel(comparePhiAngles)
        if referencePhiAngles(refIndex) == comparePhiAngles(compareIndex)
           coincidentIndex(end+1) = compareIndex; %#ok<AGROW>
        end
    end
end

compareData.correlationFunction0W0Saver = ...
    compareData.correlationFunction0W0Saver(:,coincidentIndex,:);
compareData.correlationFunction1W0Saver = ...
    compareData.correlationFunction1W0Saver(:,coincidentIndex,:);
compareData.correlationFunction2W0Saver = ...
    compareData.correlationFunction2W0Saver(:,coincidentIndex,:);

compareData.scaledR1Rates = compareData.scaledR1Rates(:,coincidentIndex,:);
compareData.r1WithPerturbationTheory = ...
    compareData.r1WithPerturbationTheory(:,coincidentIndex,:);

compareData.positionAngles = ...
    compareData.positionAngles(:,coincidentIndex,:);

compareData.configuration.('AAXXX_WARNING') = 'PHI CHANGED';
newSavingString = strrep(compare_path,'Results_orientationDependency', ...
    'Results_orientationDependency_phiCHANGED');
save(newSavingString,'-struct','compareData');
end