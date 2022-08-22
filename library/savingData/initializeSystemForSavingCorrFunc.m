function savingPath = initializeSystemForSavingCorrFunc( ...
    startDateOfSimulation,whichLipid)

savingPath = sprintf(['C:\\Users\\maxoe\\Google Drive\\Promotion' ...
    '\\Zwischenergebnisse\\relaxationRatesR1\\' ...
    '%s_%s_correlationFunctions\\'],whichLipid,startDateOfSimulation);
if ~exist(savingPath,'dir')
    mkdir(savingPath);
end

end