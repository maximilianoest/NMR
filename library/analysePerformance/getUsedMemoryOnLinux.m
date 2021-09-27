function memoryUsedInMB =  getUsedMemoryOnLinux()

% undocumented feature; see http://undocumentedmatlab.com/blog/undocumented-feature-function 
pid = feature('getpid');
[~,cmdout] = unix(['cat /proc/' num2str(pid) '/status | grep VmRSS']); % Resident Set size
 memoryUsedInMB = str2double(regexp(cmdout, '[0-9]*', 'match'))./1e6;
 
 


% more details about processes: man proc
% [~,cmdout] = unix(['cat /proc/' num2str(pid) '/status | grep VmPeak']); % Peak Virtual Memory Size
% virtualMemPeak = str2double(regexp(cmdout, '[0-9]*', 'match'))./1e3;
% 
% [~,cmdout] = unix(['cat /proc/' num2str(pid) '/status | grep VmSize']); % Virtual Memory Size
% virtualMemSize = str2double(regexp(cmdout, '[0-9]*', 'match'))./1e3;
% 
% [~,cmdout] = unix(['cat /proc/' num2str(pid) '/status | grep VmLck']); % Locked Memory Size
% virtualLckSize = str2double(regexp(cmdout, '[0-9]*', 'match'))./1e3;
% 
% [~,cmdout] = unix(['cat /proc/' num2str(pid) '/status | grep VmHWM']); % Peak Resident set size (high water mark)
% virtualHWMSize = str2double(regexp(cmdout, '[0-9]*', 'match'))./1e3;
% 
% [~,cmdout] = unix(['cat /proc/' num2str(pid) '/status | grep VmRSS']); % Resident Set size
% virtualRSSSize = str2double(regexp(cmdout, '[0-9]*', 'match'))./1e3;
% 
% 
% fprintf('MATLAB\n');
% fprintf('  Resident in RAM:   %6.0f MB\n',virtualRSSSize)
% fprintf('  Peak resident:     %6.0f MB\n',virtualHWMSize)
% fprintf('  Allocated virtual: %6.0f MB\n',virtualMemSize)
% fprintf('  Peak all. virtual: %6.0f MB\n',virtualMemPeak)
% fprintf('  Locked memory:     %6.0f MB\n',virtualLckSize)
% fprintf('\n');
end