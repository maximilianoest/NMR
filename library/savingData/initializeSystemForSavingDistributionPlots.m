function savingPath = initializeSystemForSavingDistributionPlots()

savingPath = sprintf(['C:\\Users\\maxoe\\Google Drive\\Promotion' ...
    '\\Zwischenergebnisse\\distributions\\']);
if ~exist(savingPath,'dir')
    mkdir(savingPath);
end

end
