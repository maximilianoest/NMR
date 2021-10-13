function plotCorrelationFunctionsForDifferentOrAndPos( ...
    correlationFunction,omega,dt,orientationAngles,positionAngles)

colors = getStandardizedFigure();

orientationsCount = size(orientationAngles,2);
positionsCount = size(positionAngles,2);

timeStepsCount = size(correlationFunction,3);
timeAxis = linspace(0,timeStepsCount*dt,timeStepsCount);

colorCounter = 1;
legendEntries = {};
for orientationNumber = 1:orientationsCount
    for positionNumber = 1:positionsCount
        plot(timeAxis,abs(squeeze(correlationFunction(orientationNumber ...
            ,positionNumber,:))),'Color',colors(colorCounter,:))
        colorCounter = colorCounter + 1;
        legendEntries{end+1} = ['$\theta$: ' ...
            ,num2str(orientationAngles(orientationNumber)) ...
            ,' $\varphi$: ' ...
            ,num2str(positionAngles(positionNumber))]; %#ok<AGROW>
    end
end
hold off
title(['Correlation functions q = ' num2str(omega)])
xlabel('Correlation time $\tau$ [s]')
ylabel(['$g_{' num2str(omega) '}$ ($\tau$)']);
legend(legendEntries)

end
