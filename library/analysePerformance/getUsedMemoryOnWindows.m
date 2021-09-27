function usedMemoryInMB = getUsedMemoryOnWindows()

userMemory = memory;
usedMemoryInMB = userMemory.MemUsedMATLAB/1e6; 

end

