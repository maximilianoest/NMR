function [nearestNeighboursX,nearestNeighboursY,nearestNeighboursZ ...
    ,nearestNeighbourDistancesPow3] = findNearestNeighbours( ...
    numberOfNearestNeighbours,atomIndex,relativeXPositions...
    ,relativeYPositions,relativeZPositions)
numberOfNearestNeighbours = numberOfNearestNeighbours + 1;

distances = sqrt(relativeXPositions.^2+relativeYPositions.^2 ...
    +relativeZPositions.^2);
inverseMeanDistances = 1./mean(distances,2);
inverseDistances = 0;
sumIndex = 1;

for i=1:numberOfNearestNeighbours
  smallestDistance = max(inverseMeanDistances);
  closestNeighbourId = inverseMeanDistances == smallestDistance;
  numberOfClosestElements = sum(closestNeighbourId);
  inverseDistances(sumIndex:sumIndex+numberOfClosestElements-1) = ...
      smallestDistance;
  nearestNeighbourIndex(sumIndex:sumIndex+numberOfClosestElements-1) = ...
      find(inverseMeanDistances == smallestDistance); 
  sumIndex = sumIndex + numberOfClosestElements; 
  inverseMeanDistances(closestNeighbourId) = min(inverseMeanDistances)-1;  
end

neighbourIds = false(size(numberOfNearestNeighbours));
neighbourIds(nearestNeighbourIndex) = true;
neighbourIds(atomIndex) = false;
nearestNeighbourIndex(1) = [];

nearestNeighboursX = relativeXPositions(nearestNeighbourIndex,:);
nearestNeighboursY = relativeYPositions(nearestNeighbourIndex,:);
nearestNeighboursZ = relativeZPositions(nearestNeighbourIndex,:);
nearestNeighbourDistancesPow3 = distances(nearestNeighbourIndex,:).^3;  

end
