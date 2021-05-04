clc
clearvars

%% configure the System
configuration = readConfigurationFile('config.conf');
addpath(genpath(configuration.libraryPath));

effectiveInteractionMyelinFraction = ...
    configuration.effectiveInteractionMyelinFraction;
optimizedParametersFromPaperCase1 = ...
    [0.15 0.5 2.51 1.2 14.9/effectiveInteractionMyelinFraction 16.3 1];
optimizedParametersFromPaperCase2 = ...
    [0.145 0.985 3.04 1.1 10/effectiveInteractionMyelinFraction 12.6 ...
    0.966];
optimizedParametersFromPaperCase3 = ...
    [0.14 0.59 2.57 0.86 14.7 16.0 0.99];
optimizedParametersFromPaperCase4 = ...
    [0.136 1.14 2.73 0.84 11.2 13.3 0.83];

%% calculate the real R1 Rates
