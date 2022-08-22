% This function takes the random atom sequence of a simulation given in the
% config file. Then this sequence manipulated so that it continues at the
% position where the old simulation stopped. The first part of the sequence
% is then appended to the end of the simulation. The reason for that is to
% avoid an overlap between the validation an the calculation data set.
% IN CONFIG FILE: (path2OldSequence)
% EXAMPLE:
% sequence of old simulation = [ 1 2 3 4 5 6 7 8 9 10]
% in this old simulation 5 atoms where calculated. Thus, the new sequence
% is:
% [6 7 8 9 10 1 2 3 4 5]

function sequenceForNonOverlappingDatasets = ...
    reloadOldAtomSequenceToAvoidOverlappingDataSets()
configuration = readConfigurationFile('config.txt');

if ~isfield(configuration,'path2OldSequence')
    error('reloadOldSequenceError:configurationNotFound', ...
        ['The configuration on "path2OldSequence" cannot be found.' ...
        ' Check the config file.'])
end

if exist(configuration.path2OldSequence,'file')
    results = load(configuration.path2OldSequence);
    atomSequence = results.randomSequenceOfAtoms;
    atomCounter = results.atomCounter;
    sequenceForNonOverlappingDatasets = [atomSequence(atomCounter+1:end) ...
        atomSequence(1:atomCounter)];
    
    fprintf('    The old atom sequence starts with\n    ');
    fprintf('%i ',atomSequence(1:10));
    fprintf('\n');
    fprintf('    The atom counter were %i\n',atomCounter);
    fprintf('    The old sequence at the atom counter position is:\n    ');
    fprintf('%i ',atomSequence(atomCounter+1:atomCounter+11));
    fprintf('\n');
    fprintf('    Thus, the sequence of this new simulation starts with:\n    ');
    fprintf('%i ',sequenceForNonOverlappingDatasets(1:10));
    fprintf('\n');
    
else
    error('reloadOldSequenceError:fileNotFound', ...
        'The file %s cannot be found',configuration.path2OldSequence);
end


end

