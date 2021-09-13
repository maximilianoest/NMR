function sizeInGigaByte =  getSize(variable)  %#ok<INUSD>
   information = whos('variable');
   mem = memory;
   mem.MemUsedMATLAB/1e9
end