function [atomIndices] = getIndicesOfSnippet(moleculeSize,headSize ...
    ,offset,snippetSize,numberOfMolecules)

atomIndices = zeros(numberOfMolecules,snippetSize);

for moleculeNr = 1:numberOfMolecules
        startingPoint = headSize+offset+(moleculeNr-1)*moleculeSize;
        endPoint = startingPoint+snippetSize;
        atomIndices(moleculeNr,:) = startingPoint:endPoint-1;
end
end