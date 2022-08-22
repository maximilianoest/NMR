function showResultsOfScaledUpR1(scaledR1_path)

initializeFigure();
legendEntries =[];
data = load(scaledR1_path);
scaledUpR1_theta_phi = data.scaledUpR1_theta_phi;
orientationAngles = data.orientationAngles;

% theta and phi dependent R1
for posNr = 1:length(data.positionAngles)
    plot(orientationAngles,scaledUpR1_theta_phi(:,posNr));
    legendEntries{end+1} = sprintf('$\\varphi$: %.2f', ...
        data.positionAngles(posNr)); %#ok<AGROW>
end

% theta dependent R1 (= effective theta dependent R1)
effectiveThetaDependentR1 = squeeze(mean(scaledUpR1_theta_phi,2))'

plot(orientationAngles,effectiveThetaDependentR1,'--');
legendEntries{end+1} = ('effective $\theta$-dependent');

overallR1 = mean(effectiveThetaDependentR1) %#ok<NASGU>
if effectiveThetaDependentR1(1) > effectiveThetaDependentR1(end)
    location = 'southwest';
else
    location = 'northwest';
end


legend(legendEntries,'location',location);
title(sprintf('Orientation-dependent R$_1$ (%s, Matl.: %s, Grom.: %s)', ...
    data.whichLipid,data.matlabSimulationDate,data.gromacsSimulationDate));
xlabel('$\theta$ (angle of axon to B$_0$)');
ylabel('R$_1$ [Hz]');

savingPath = initializeSystemForSavingR1Plots();
fieldstrengthString = strrep(num2str(data.fieldStrength),'.','');
savingName = sprintf('%s_r1%s_%sTesla_AngleDependentUpR1', ...
    data.whichLipid,data.matlabSimulationDate,fieldstrengthString);
print(gcf,[savingPath savingName],'-dpng','-r300');

end

