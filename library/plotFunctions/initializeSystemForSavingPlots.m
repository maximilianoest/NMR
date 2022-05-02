function savingPath = initializeSystemForSavingPlots(kindOfSaving ...
    ,whichLipid)

savingPath = sprintf(['C:\\Users\\maxoe\\Google Drive\\Promotion' ...
    '\\Zwischenergebnisse\\%s_%s\\%s_%s\\'] ...
    ,kindOfSaving,whichLipid,datestr(date,'yyyymmdd'),whichLipid);
if ~exist(savingPath,'dir')
    mkdir(savingPath);
end

end
