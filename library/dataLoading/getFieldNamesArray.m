function [fieldNamesArray] = getFieldNamesArray(fieldNamesToLoad)

fieldNamesCell = split(fieldNamesToLoad,",");

for fieldNameNumber = 1:length(fieldNamesCell)
    fieldName = string(fieldNamesCell{fieldNameNumber});
    fieldNamesArray(fieldNameNumber) = fieldName; %#ok<AGROW>
end

if isempty(fieldNamesArray)
   warning('No Field Names where given'); 
end

end