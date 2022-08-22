%#ok<*AGROW>
function [results_path] = calculateR1ForNNPostProcessing( ...
    resultsCorrFunc_path,fieldStrength,OFFSETSUPPRESSIONREGION, ...
    CUTOFFFRACTION,PLOTCORRELATIONFUNCTIONS,savingAdditional)

% constants stuff
constants = readConstantsFile("../../constants.txt");
dipolDipolConstant = 3/4*(constants.vaccumPermeability/(4*pi) ...
    *constants.hbar*constants.gyromagneticRatioOfHydrogenAtom^2)^2 ...
    /(constants.nanoMeter^6);
gyromagneticRatio = constants.gyromagneticRatioOfHydrogenAtom;

% contribution of magnetic field strength
omega0 = gyromagneticRatio*fieldStrength;

% results from NN dependent correlation determination
results = load(resultsCorrFunc_path);
whichLipid = results.whichLipid;
atomCounter = results.atomCounter;
deltaT = results.samplingFrequency;
correlationFunctions1W0 = double(results.correlationFunction1W0Saver);
correlationFunctions2W0 = double(results.correlationFunction2W0Saver);
corrFuncLength = size(correlationFunctions1W0,4);
gromacsFileName = results.fileName;

matlabSimulationDate = results.startDateOfSimulation;
gromacsSimulationDate = results.simulationDate;
nearestNeighbourCases = results.nearestNeighbourCases;
orientationAngles = rad2deg(results.orientationAngles);
positionAngles = rad2deg(results.positionAngles);

offsetSuppressionRegion = round(OFFSETSUPPRESSIONREGION(1) ...
    *corrFuncLength):round(OFFSETSUPPRESSIONREGION(2)*corrFuncLength);
remainingRegion = [1:round(corrFuncLength*CUTOFFFRACTION)-1];
for nnCaseNr = 1:length(nearestNeighbourCases)
    for orientationNr = 1:length(orientationAngles)
        for positionNr = 1:length(positionAngles)
            
            % offset suppression due to static part of Hamiltonian
            correlationFunctions1W0(nnCaseNr,orientationNr,positionNr,:) = ...
                suppressOffsetOfCorrFunc(correlationFunctions1W0( ...
                nnCaseNr,orientationNr,positionNr,:), ...
                offsetSuppressionRegion);
            correlationFunctions2W0(nnCaseNr,orientationNr,positionNr,:) = ...
                suppressOffsetOfCorrFunc(correlationFunctions2W0( ...
                nnCaseNr,orientationNr,positionNr,:), ...
                offsetSuppressionRegion);
            
            % cut off part of correlation function to avaoid rise of
            % correlation function with increasing tau caused by FT to
            % determine correlation function.
            cuttedCorrelationFunctions1W0(nnCaseNr,orientationNr, ...
                positionNr,:) = correlationFunctions1W0(nnCaseNr, ...
                orientationNr,positionNr,remainingRegion);
            cuttedCorrelationFunctions2W0(nnCaseNr,orientationNr, ...
                positionNr,:) = correlationFunctions2W0(nnCaseNr, ...
                orientationNr,positionNr,remainingRegion);
        end
    end
end

% calculate R1 from cutted and offset suppressed correlation function
lags = size(cuttedCorrelationFunctions1W0,4);
for nnCaseNr = 1:length(nearestNeighbourCases)
    for orientationNr = 1:length(orientationAngles)
        for positionNr = 1:length(positionAngles)
            r1_theta_phi_NN(orientationNr,positionNr, ...
                nnCaseNr) = calculateR1FromCorrFuncWithExplIntegration( ...
                squeeze(cuttedCorrelationFunctions1W0(nnCaseNr,orientationNr, ...
                positionNr,:))',squeeze(cuttedCorrelationFunctions2W0(nnCaseNr, ...
                orientationNr,positionNr,:))',deltaT,omega0, ...
                lags,dipolDipolConstant);
        end
    end
end

savingPath = initializeSystemForSavingR1();
fieldstrengthString = strrep(num2str(fieldStrength),'.','');
savingName = sprintf('%s_%s_%sTesla_relaxationRates%s', ...
    whichLipid,matlabSimulationDate,fieldstrengthString,savingAdditional);
results_path = [savingPath savingName '.mat'];
save(results_path,'r1_theta_phi_NN','whichLipid','fieldStrength', ...
    'matlabSimulationDate','gromacsSimulationDate','orientationAngles','positionAngles', ...
    'nearestNeighbourCases','atomCounter','gromacsFileName');

if PLOTCORRELATIONFUNCTIONS
    timeAxis = [0:deltaT:(lags-1)*deltaT];
    for nnCaseNr = 1:length(nearestNeighbourCases)
        initializeFigure();
        legendEntries = {};
        for orientationNr = 1:length(orientationAngles)
            for positionNr = 1:length(positionAngles)
                plot(timeAxis,squeeze(cuttedCorrelationFunctions1W0(nnCaseNr, ...
                    orientationNr,positionNr,:) ...
                    /cuttedCorrelationFunctions1W0(nnCaseNr, ...
                    orientationNr,positionNr,1)));
                legendEntries{end+1} = sprintf(['first order, $\\theta$: %.2f ' ...
                    '$\\varphi$: %.2f'], ...
                    orientationAngles(orientationNr), ...
                    positionAngles(positionNr));
                plot(timeAxis,squeeze(cuttedCorrelationFunctions2W0(nnCaseNr, ...
                    orientationNr,positionNr,:) ...
                    /cuttedCorrelationFunctions2W0(1, ...
                    orientationNr,positionNr,1)));
                legendEntries{end+1} = sprintf(['second order $\\theta$: %.2f ' ...
                    '$\\varphi$: %.2f'], ...
                    orientationAngles(orientationNr), ...
                    positionAngles(positionNr));
                
            end
        end
        title(sprintf(['Cleaned correlation functions (NN:%i %s)'], ...
            nearestNeighbourCases(nnCaseNr),results.whichLipid));
        xlabel('Correlation time [s]');
        ylabel('Normalized amplitude [a.u.]')
        legend(legendEntries);
        savingPathCorrFunc = initializeSystemForSavingCorrFunc( ...
            matlabSimulationDate,whichLipid);
        savingName = sprintf('%s_normOffsetSupprAndCuttedCorrFunc_NN%i', ...
            results.startDateOfSimulation, ...
            nearestNeighbourCases(nnCaseNr));
        print(gcf,[savingPathCorrFunc savingName],'-dpng','-r300');
    end
end

end