function offsetSuppressedCorrFunc = suppressOffsetOfCorrFunc( ...
    correlationFunction,offsetSuppressionRegion)
% This function suppresses the offset of the corrlelation function for
% given region where the offset should be nearly constant.

    offsetSuppressedCorrFunc = correlationFunction...
                - mean(correlationFunction(offsetSuppressionRegion));
end
