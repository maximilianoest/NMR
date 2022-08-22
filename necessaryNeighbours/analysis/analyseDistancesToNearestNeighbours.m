close all

results = load("C:\Users\maxoe\Google Drive\Promotion\Results\nearestNeighboursAnalysis\Server\20211013_Results_relevantNearestNeighbours_water_H_50ns_05ps_wh.mat");
meanPositions = results.meanPositions;
edges = [0:0.01:1];

for atom = 1:100
    randomAtom = randi([500,size(meanPositions,1)-4000],1);
    relativePositions = meanPositions - meanPositions(randomAtom,:);
    relativeDistances = sqrt(sum(relativePositions.^2,2));
    inverseMeanDistances = 1./relativeDistances;
    inverseDistances = 0;
    sumIndex = 1;
    numberOfNearestNeighbours = 3000;

    for nearestNeighbour=1:numberOfNearestNeighbours
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
    inverseDistances(1) = [];
    distances = 1./inverseDistances;
    for counter = 1:length(edges)-1
        barCount(atom,counter) = sum(distances < edges(counter+1) ...
            & distances >= edges(counter));
    end
    
%     figure(1)
%     quantity = histogram(distances,edges);
%     quantity = quantity.Data;
%     grid minor
%     pause(2)
end

% figure(1)
% xlabel('Distance [nm]')
% xlabel('Frequency')
% xlabel('Distance [nm]')
% ylabel('Frequency')
% title(['Frequency of distances for atom ' num2str(randomAtom) ])

averageNearestNeighbours = mean(barCount,1);
figure(2)
bar(edges(1:end-1),averageNearestNeighbours)
title('Number of nearest neighbours at a certain distance.')
xlabel('Distance [nm]')
ylabel('Average frequency')
grid minor


figure(3)
plot(edges,1./edges.^3,'LineWidth',1.5)
title('Influence of inverse distance on spectral density')
xlabel('Distance [nm]')
ylabel('1/Distance^3')
grid minor





% neighbourIds = false(size(numberOfNearestNeighbours));
% neighbourIds(nearestNeighbourIndex) = true;
% neighbourIds(atomIndex) = false;
% nearestNeighbourIndex(1) = [];
% nearestNeighboursX = relativeXPositions(nearestNeighbourIndex,:);
% nearestNeighboursY = relativeYPositions(nearestNeighbourIndex,:);
% nearestNeighboursZ = relativeZPositions(nearestNeighbourIndex,:);
% nearestNeighbourDistancesPow3 = distances(nearestNeighbourIndex,:).^3;  