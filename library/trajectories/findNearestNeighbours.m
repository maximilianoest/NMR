function [nearestNeighboursX,nearestNeighboursY,nearestNeighboursZ ...
    ,nearestNeighbourDistancesPow3] = findNearestNeighbours( ...
    numberOfNearstNeighbours,atomIndex,relativeXPositions...
    ,relativeYPositions,relativeZPositions)
% this function searches for the closest nearest neighbours and
% calculates the distance in X,Y and Z direction.

distances = sqrt(relativeXPositions.^2+relativeYPositions.^2 ...
    +relativeZPositions.^2);
inverseMeanDistances = 1./mean(distances,2);
inverseDistances = 0;
sumIndex = 1;

for i=1:numberOfNearstNeighbours
  smallestDistance = max(inverseMeanDistances);
  closestNeighbourId = inverseMeanDistances == smallestDistance;
  numberOfClosestElements = sum(closestNeighbourId);
  inverseDistances(sumIndex:sumIndex+numberOfClosestElements-1) = ...
      smallestDistance;
  nearestNeighbourIndex(sumIndex:sumIndex+numberOfClosestElements-1) = ...
      find(inverseMeanDistances == smallestDistance); 
  sumIndex = sumIndex + numberOfClosestElements; 
  inverseMeanDistances(closestNeighbourId) = min(inverseMeanDistances);  
end

neighbourIds = false(size(distances));
neighbourIds(nearestNeighbourIndex) = true;
neighbourIds(atomIndex) = false;

nearestNeighboursX = relativeXPositions(neighbourIds,:);
nearestNeighboursY = relativeYPositions(neighbourIds,:);
nearestNeighboursZ = relativeZPositions(neighbourIds,:);
nearestNeighbourDistancesPow3 = distances(neighbourIds,:).^3;  
end
