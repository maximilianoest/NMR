function [constants] = readConstantsFile(path2ConstantsFile)

constants = {};
constantFileId = fopen(path2ConstantsFile);
data = textscan(constantFileId, '%s %s','Delimiter','=');
fclose(constantFileId);

for element = 1:size(data{1},1)
    constantName = data{1}{element};
    try
        constantData = data{2}{element};
        if isnan(str2double(constantData))
            warning(['The constant "' constantName '" can not be read.' ...
                'Please check if you gave a numerical value or a string.'])
        else 
            constants.(constantName) = str2double(data{2}{element});
        end
        
    catch 
        warning(['For the constant variable "' constantName ...
            '" no data are given. Take a look at the configuration file.']);
    end
end


end
