function [trimmedData] = trimmR1DataWithPercentile(calculatedR1Rates ...
    ,lowerPercentile,upperPercentile,orientationCount,positionCount)

for orientationNumber = 1:orientationCount
    for positionNumber = 1:positionCount
        data = squeeze(calculatedR1Rates(orientationNumber ...
            ,positionNumber,:))';
        trimmedIndices = data > prctile(data,lowerPercentile) ...
            & data < prctile(data,upperPercentile);
        trimmedData(orientationNumber,positionNumber,:) ...
            = data(trimmedIndices);  %#ok<AGROW>
    end
end

end
