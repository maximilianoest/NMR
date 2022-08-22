function savingPath = initializeSystemForSavingR1Plots()

savingPath = sprintf(['C:\\Users\\maxoe\\Google Drive\\Promotion' ...
    '\\Zwischenergebnisse\\relaxationRatesR1\\Plots\\']);
if ~exist(savingPath,'dir')
    mkdir(savingPath);
end

end