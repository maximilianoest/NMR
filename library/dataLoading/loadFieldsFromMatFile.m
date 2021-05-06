function [dataFields] = loadFieldsFromMatFile(path2Data ...
    ,dataFieldNameStrings)

dataFields = struct();
data = load(path2Data);

for fieldName = dataFieldNameStrings
    try
        dataFields.(fieldName) = data.(fieldName);
    catch
        warning(['Field name "' char(fieldName) '" does not exist in' ...
            ' the data you try to read. This data field is skipped.']);
    end
end

end
