function [correlationFunctions1W0,correlationFunctions2W0] = ...
    getCorrelationFunctionsFromSimulationData(data)

correlationFunctions1W0 = data.correlationFunction1W0Saver;
correlationFunctions2W0 = data.correlationFunction2W0Saver;
atomCount = data.atomCounter;
correlationFunctions1W0 = ...
    squeeze(correlationFunctions1W0(:,:,1:atomCount,:));
correlationFunctions2W0 = ...
    squeeze(correlationFunctions2W0(:,:,1:atomCount,:));

end
