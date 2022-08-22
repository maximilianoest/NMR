function compareNNDependentScaledR1(r1ScaledUp_path,saving)
load(r1ScaledUp_path); %#ok<LOAD>

% 1. plot theta-dependent and overall R1 values
initializeFigure();
legendEntries = {};
for oriNr = 1:size(r1ScaledUp_theta_phi_NN,1)
    for posNr = 1:size(r1ScaledUp_theta_phi_NN,2)
        plot(nearestNeighbourCases, ...
            squeeze(r1ScaledUp_theta_phi_NN(oriNr,posNr,:)),'+-');
        legendEntries{end+1} = sprintf('R1, $\\theta$: %4.2f $\\varphi$: %4.2f', ...
            orientationAngles(oriNr),positionAngles(posNr));
    end
end
for oriNr = 1:size(r1ScaledUp_theta_NN,1)
    plot(nearestNeighbourCases, ...
        squeeze(r1ScaledUp_theta_NN(oriNr,:)),'.-');
    legendEntries{end+1} = sprintf('Effective R1, $\\theta$: %4.2f', ...
        orientationAngles(oriNr));
end
plot(nearestNeighbourCases,r1ScaledUp_NN,'*-');
legendEntries{end+1} = sprintf('Overall R1');

legend(legendEntries,'location','east');
xlabel('Nearest neighbours')
ylabel('R$_1$ [Hz]')
title(sprintf('Scaled up R$_1$ (%s)',whichLipid));
legend(legendEntries)

if saving
    savingPath = initializeSystemForSavingR1Plots();
    savingName = sprintf('%s_r1%s_scalFac%s_scaledUpR1_NNDependent_%s', ...
        whichLipid,r1MatlabSimulationDate,scalFacMatlabSimulationDate);
    print(gcf,[savingPath savingName],'-dpng','-r300');
end

% 2. plot theta-dependent R1 shift
initializeFigure();
legendEntries={};

for oriNr = 1:size(r1ScaledUp_theta_NN,1)
    plot(nearestNeighbourCases, ...
        squeeze(r1ScaledUp_theta_NN(oriNr,:)- r1ScaledUp_theta_NN(1,:)));
    legendEntries{end+1} = sprintf('$\\theta$: %4.2f', ...
        orientationAngles(oriNr));
end
legend(legendEntries,'location','east')
xlabel('Nearest neighbours')
ylabel('R$_1$ shift [Hz]')
title(sprintf('$\\theta$-dependent scaled up R$_1$ shift (%s)' ...
    ,whichLipid));
legend(legendEntries,'Location','East')

if saving
    savingName = sprintf('%s_r1%s_scalFac%s_scaledUpR1Shift_thetaNNDependent_%s', ...
        whichLipid,r1MatlabSimulationDate,scalFacMatlabSimulationDate);
    print(gcf,[savingPath savingName],'-dpng','-r300');
end


% 3. plot difference in overall R1 to highest NN case
% overallDifference = overallPredictedR1 - overallPredictedR1(:,:,1);
% 
% initializeFigure();
% legendEntries = {};
% plot(nearestNeighbourCases,squeeze(overallDifference(refOriNr,refPosNr,:)),'*-');
%         legendEntries{end+1} = sprintf('Ref.: $\\theta$: %.2f, $\\varphi$: %.2f',orientationAngles(refOriNr),positionAngles(refPosNr));
% 
% for refOriNr = 1:orientationsCount
%     for refPosNr = 1:positionsCount
%         plot(nearestNeighbourCases,squeeze(overallDifference(refOriNr,refPosNr,:)),'*-');
%         legendEntries{end+1} = sprintf('Ref.: $\\theta$: %.2f, $\\varphi$: %.2f',orientationAngles(refOriNr),positionAngles(refPosNr));
%     end
% end
% legend(legendEntries,'location','east')
% 
% xlabel('Nearest neighbours')
% ylabel('Difference [Hz]')
% title(sprintf('Difference in predicted overall R$_1$ to higher NN (%s)',whichLipid));
% 
% if saving
%     savingName = sprintf('%s_predictedR1DiffrenceToHigh_NNDependent_%s' ...
%         ,results.startDateOfSimulation,whichLipid);
%     print(gcf,[savingPath savingName],'-dpng','-r300');
% end




end

% r1ScaledUp_path='C:\Users\maxoe\Google Drive\Promotion\Zwischenergebnisse\relaxationRatesR1\PLPC_r120220425_scalFac20220425_3Tesla_relaxationRatesScaledUp.mat'