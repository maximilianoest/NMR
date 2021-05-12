function [trimmedData] = trimmR1DataWithPercentile(calculatedR1Rates ...
    ,percentile)

percentile = prctile(calculatedR1Rates,percentile,[1,2]);

disp('test');


end
