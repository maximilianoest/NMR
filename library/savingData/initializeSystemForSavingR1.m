function savingPath = initializeSystemForSavingR1()

savingPath = sprintf(['C:\\Users\\maxoe\\Google Drive\\Promotion' ...
    '\\Zwischenergebnisse\\relaxationRatesR1\\']);
if ~exist(savingPath,'dir')
    mkdir(savingPath);
end

end