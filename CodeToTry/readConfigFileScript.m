clc
clear all

configuration = {};
configFileId = fopen('config.conf');
data = textscan(configFileId, '%s %s','Delimiter','=');
fclose(configFileId);

for element = 1:size(data{1},1)
    configurationName = data{1}{element};
    try
        configuration.(configurationName) = str2double(data{2}{element});
    catch 
        warning(['For the configuration variable "' configurationName ...
            '" no data are given. Take a look at the configuration file.']);
    end
end

