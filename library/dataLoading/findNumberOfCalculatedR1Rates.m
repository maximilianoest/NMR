function [atomCount] = findNumberOfCalculatedR1Rates(r1Array)

atomsWithCalculatedR1Rates = (r1Array > 0.0001);
atomsWithCalculatedR1Rates = logical(squeeze(mean(mean( ...
    atomsWithCalculatedR1Rates,1),2))');
atomCount = sum(atomsWithCalculatedR1Rates);

end
