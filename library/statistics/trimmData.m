function [trimmedData] = trimmData(data ...
    ,lowerPercentile,upperPercentile)

trimmedIndices = data > prctile(data,lowerPercentile) ...
    & data < prcentile(data,upperPercentile);
trimmedData = data(trimmedIndices);


end
