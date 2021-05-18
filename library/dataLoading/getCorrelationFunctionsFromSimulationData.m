function [correlationFunctions1W0,correlationFunctions2W0] = ...
    getCorrelationFunctionsFromSimulationData(data)

atomCount = data.atomCounter;
correlationFunctions1W0 = data.correlationFunction1W0Saver;
correlationFunctions2W0 = data.correlationFunction2W0Saver;
atomIndex = data.atomIndex;
atomIndex = atomIndex(1:atomCount);
correlationFunctions1W0 = ...
    squeeze(correlationFunctions1W0(:,:,atomIndex,:));
correlationFunctions2W0 = ...
    squeeze(correlationFunctions2W0(:,:,atomIndex,:));

end
