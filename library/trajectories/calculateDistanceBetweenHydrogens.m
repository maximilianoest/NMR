function [distance] = ...
    calculateDistanceBetweenHydrogens(trajectory)
firstHydrogenInWaterTrajectory = trajectory(1:2:end,:);
secondHydrogenInWaterTrajectory = trajectory(2:2:end,:);

if size(firstHydrogenInWaterTrajectory,1) ~= ...
        size(secondHydrogenInWaterTrajectory,1)
    firstHydrogenInWaterTrajectory = ...
        firstHydrogenInWaterTrajectory(1:end-1,:);
end
distance = firstHydrogenInWaterTrajectory-secondHydrogenInWaterTrajectory;

end
