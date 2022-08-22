function savingPath = initializeSystemForSavingScalingFactors()

savingPath = sprintf(['C:\\Users\\maxoe\\Google Drive\\Promotion' ...
    '\\Zwischenergebnisse\\scalingFactors\\']);
if ~exist(savingPath,'dir')
    mkdir(savingPath);
end

end
