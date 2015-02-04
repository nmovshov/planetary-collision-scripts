function fnl = load_fnl(filename)
% Load node list data from file and parse out to a table.
%
% This is a slimmed-down version of the method in ahelpers.py. This function
% loads all node lists in the file into a single table (v >= R2014b) or array.

% Unzip if needed
[~,~,ext] = fileparts(filename);
if strcmp(ext,'.gz')
    filename = gunzip(filename,tempdir);
    filename = filename{1};
    xxxxxxxx = onCleanup(@()delete(filename)); % leave no file behind...
end

% Load and distribute data
raw = importdata(filename,' ',headcount(filename));
varNames = {'id', 'eos_id', 'x', 'y', 'z', 'vx', 'vy', 'vz', 'm', 'rho',...
            'P', 'T', 'U', 'hmin', 'hmax'};
if verLessThan('matlab','8.4.0')
    fnl = raw.data;
else
    fnl = array2table(raw.data, 'VariableNames', varNames);
end
