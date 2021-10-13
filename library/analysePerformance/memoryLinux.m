function memoryLinux()
% Brief memory usage overview

[~,cmdout] = unix('free | grep Mem');

stats = str2double(regexp(cmdout, '[0-9]*', 'match'));
stats = stats./1e3; % kB to MB

totalMemory    = stats(1);
usedMemory     = stats(2);
freeMemory     = stats(3);
sharedMemory   = stats(4);
bufferedMemory = stats(5);
cashedMemory   = stats(6);
freeAndCashMemory = freeMemory+cashedMemory;

% undocumented feature; see http://undocumentedmatlab.com/blog/undocumented-feature-function 
pid = feature('getpid');

% more details about processes: man proc
[~,cmdout] = unix(['cat /proc/' num2str(pid) '/status | grep VmPeak']); % Peak Virtual Memory Size
virtualMemPeak = str2double(regexp(cmdout, '[0-9]*', 'match'))./1e3;

[~,cmdout] = unix(['cat /proc/' num2str(pid) '/status | grep VmSize']); % Virtual Memory Size
virtualMemSize = str2double(regexp(cmdout, '[0-9]*', 'match'))./1e3;

[~,cmdout] = unix(['cat /proc/' num2str(pid) '/status | grep VmLck']); % Locked Memory Size
virtualLckSize = str2double(regexp(cmdout, '[0-9]*', 'match'))./1e3;

[~,cmdout] = unix(['cat /proc/' num2str(pid) '/status | grep VmHWM']); % Peak Resident set size (high water mark)
virtualHWMSize = str2double(regexp(cmdout, '[0-9]*', 'match'))./1e3;

[~,cmdout] = unix(['cat /proc/' num2str(pid) '/status | grep VmRSS']); % Resident Set size
virtualRSSSize = str2double(regexp(cmdout, '[0-9]*', 'match'))./1e3;

% Formated output 
fprintf('\n');
fprintf('SYSTEM\n');
fprintf('  Total memory:  %6.0f MB\n',totalMemory)
fprintf('  Used memory:   %6.0f MB\n',usedMemory)
fprintf('  Free memory:   %6.0f MB\n',freeMemory)
fprintf('  Shared memory: %6.0f MB\n',sharedMemory)
fprintf('  Buffer memory: %6.0f MB\n',bufferedMemory)
fprintf('  Cashed memory: %6.0f MB\n',cashedMemory)
fprintf('  Free + cashed: %6.0f MB\n',freeAndCashMemory)
fprintf('MATLAB\n');
fprintf('  Resident in RAM:   %6.0f MB\n',virtualRSSSize)
fprintf('  Peak resident:     %6.0f MB\n',virtualHWMSize)
fprintf('  Allocated virtual: %6.0f MB\n',virtualMemSize)
fprintf('  Peak all. virtual: %6.0f MB\n',virtualMemPeak)
fprintf('  Locked memory:     %6.0f MB\n',virtualLckSize)
fprintf('\n');
end