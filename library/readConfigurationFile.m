function [configuration] = readConfigurationFile(path2ConfigFile)

configuration = {};
configFileId = fopen(path2ConfigFile);
data = textscan(configFileId, '%s %s','Delimiter','=');
fclose(configFileId);

for element = 1:size(data{1},1)
    configurationName = data{1}{element};
    try
        configurationInformation = data{2}{element};
        if isnan(str2double(configurationInformation))
            configuration.(configurationName) = configurationInformation;
        else
            configuration.(configurationName) = ...
                str2double(configurationInformation);
        end 
        
    catch 
        warning(['For the configuration variable "' configurationName ...
            '" no data are given. Take a look at the configuration file.']);
    end
end

end
