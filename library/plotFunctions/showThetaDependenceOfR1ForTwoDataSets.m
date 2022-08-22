function showThetaDependenceOfR1ForTwoDataSets(firstDataSet_path, ...
    secondDataSet_path,saving)

firstSet = load(firstDataSet_path);
secondSet = load(secondDataSet_path);


firstR1_theta_phi = firstSet.r1_theta_phi_NN(:,:,1);
secondR1_theta_phi = secondSet.r1_theta_phi_NN(:,:,1);

firstR1_theta = squeeze(mean(firstR1_theta_phi,2));
secondR1_theta = squeeze(mean(secondR1_theta_phi,2));

orientationAngles = firstSet.orientationAngles;

initializeFigure();
legendEntries = {};
plot(orientationAngles,firstR1_theta);
legendEntries{end+1} = sprintf('%i atoms',firstSet.atomCounter);
plot(orientationAngles,secondR1_theta);
legendEntries{end+1} = sprintf('%i atoms',secondSet.atomCounter);
legend(legendEntries);
xlabel('$\theta$');
ylabel('Relaxation rate R$_1$ [Hz]');
title(sprintf('%s comparison $\\theta$-dependent R$_1$', ...
    firstSet.whichLipid));

if saving
    savingPath = initializeSystemForSavingR1Plots();
    savingName = sprintf('%s_%s_thetaDependentR1ForDifferentAtomCounts_%s' ...
        ,firstSet.matlabSimulationDate, ...
        secondSet.matlabSimulationDate,firstSet.whichLipid);
    print(gcf,[savingPath savingName],'-dpng','-r300');
end

end
