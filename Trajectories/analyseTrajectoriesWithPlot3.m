clc
clear all
close all


load("C:\Users\maxoe\Google Drive\Promotion\Data\Lipids\DOPS\20220110_DOPS_TIP4_Bilayer_50H20\20220110_DOPS_TIP4_Bilayer_50water_lipid_H_whole_dt2ps_simTime1000ns_reducedTrajectories.mat");
lastFrame = size(singleFrameTrajectoryX,1);
hydrogenAtoms = 200;
frameNr = lastFrame;
figure;
hold on
plot3(singleFrameTrajectoryX(frameNr,1:hydrogenAtoms) ...
    ,singleFrameTrajectoryY(frameNr,1:hydrogenAtoms) ...
    ,singleFrameTrajectoryZ(frameNr,1:hydrogenAtoms) ...
    ,'.','LineWidth',2);

plot3(singleFrameTrajectoryZ(frameNr,1:hydrogenAtoms) ...
    ,singleFrameTrajectoryY(frameNr,1:hydrogenAtoms) ...
    ,singleFrameTrajectoryX(frameNr,1:hydrogenAtoms) ...
    ,'.','LineWidth',2);

rotationAxis = [0 1 0];
rotationAngle = deg2rad(90);
rotationMatrix = get3DRotationMatrix(rotationAngle,rotationAxis);
[rotatedX,rotatedY,rotatedZ] = rotateTrajectoriesWithRotationMatrix( ...
    rotationMatrix,singleFrameTrajectoryZ(frameNr,1:hydrogenAtoms) ...
    ,singleFrameTrajectoryY(frameNr,1:hydrogenAtoms) ...
    ,singleFrameTrajectoryX(frameNr,1:hydrogenAtoms));
plot3(rotatedX,rotatedY,rotatedZ,'*','LineWidth',1);

hold off
axis([-10 10 -10 10 -10 10])
grid on
xlabel('X')
ylabel('Y')
zlabel('Z')
set(gca,'FontSize',12)
set(gcf,'WindowState','maximized')
legend('My frame','GROMACS frame','Rotated with function');

% figure;
% grid minor
% for frameNr = 1:lastFrame
%     plot3(singleFrameTrajectoryZ(frameNr,:) ...
%         ,singleFrameTrajectoryY(frameNr,:) ...
%         ,singleFrameTrajectoryX(frameNr,:),'.');
%     xlabel('X')
%     ylabel('Y')
%     zlabel('Z')
%     axis([-10 10 -10 10 -10 10])
%     set(gca,'FontSize',12)
%     set(gcf,'WindowState','maximized')
%     title('Original frame from GROMACS')
%     pause(0.1)
% end

    
    
    
    
    
    
    