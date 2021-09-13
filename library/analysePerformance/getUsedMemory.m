function usedMemoryInGB = getUsedMemory()

variableInformation = evalin('base','whos');
if size(variableInformation) > 0
    for elementNr = 1:size(variableInformation)
       memoryUsedByVariables(elementNr) = variableInformation( ...
           elementNr).bytes;
    end
    usedMemoryInGB = sum(memoryUsedByVariables)/(2^30);
else
    usedMemoryInGB = 0;
end

end
