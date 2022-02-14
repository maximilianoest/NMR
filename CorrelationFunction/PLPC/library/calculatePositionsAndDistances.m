function [relativeXPositions,relativeYPositions,relativeZPositions ...
    ,distances] = calculatePositionsAndDistances(trajX,trajY ...
    ,trajZ,atomNumber)
% This function calculates the relative x, y and z positions of
% the hydrogen atoms with respect to the hydrogen atom of interest.
% Additionally the distances between these atoms are calculated.

% positions of obeserved H in world coordinate system
xPosition = trajX(atomNumber,:);
yPosition = trajY(atomNumber,:);
zPosition = trajZ(atomNumber,:);

% positions of other particles relative to observed H
relativeXPositions = bsxfun(@minus,trajX,xPosition);
relativeYPositions = bsxfun(@minus,trajY,yPosition);
relativeZPositions = bsxfun(@minus,trajZ,zPosition);
distances = sqrt(relativeXPositions.^2+relativeYPositions.^2 ... 
    +relativeZPositions.^2);

    end