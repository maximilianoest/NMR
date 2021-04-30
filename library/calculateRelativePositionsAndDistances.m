function [trajX,trajY,trajZ,distances] = ...
    calculateRelativePositionsAndDistances(trajX,trajY,trajZ,atomNumber)

trajX = bsxfun(@minus,trajX,trajX(atomNumber,:));
trajY = bsxfun(@minus,trajY,trajY(atomNumber,:));
trajZ = bsxfun(@minus,trajZ,trajZ(atomNumber,:));
distances = sqrt(trajX.^2+trajY.^2+trajZ.^2);

end