% GENERATE_PSET
% Generates a point-set cell array from a zipped data file
%
% pset = generate_pset(zipfile, nsamples, ntrials)
%
% where
%
% ZIPFILE is the path name or URL of the zipped data file containing the
% circuit measurements. Default:
% http://www.cs.tut.fi/~timhome/tutorials/circuits/files/tutorial-circuits-data.zip
%
% NSAMPLES is the number of data samples to use in the analysis. NSAMPLES
% must be smaller or equal than 3000. Default: 3000.
%
% NTRIALS is the number of data trials to use. At most 182 trials can be
% used. Default: 182.
%
% PSET is a cell array of dimensions (2 x NTRIALS) containing the
% measurement time-series.
%
% See also: tutorial_circuits_analysis, tutorial_circuits_figures

% Description: Point-set generation
% Documentation: tutorial_circuits.txt
% Author: German Gomez-Herrero


function pset = generate_pset(zipfile, nsamples, ntrials)

MAX_TRIALS = 182;
MAX_SAMPLES = 3000;

if nargin < 3 || isempty(ntrials), ntrials = MAX_TRIALS; end
if nargin < 2 || isempty(nsamples), nsamples = MAX_SAMPLES; end
if nargin < 1 || isempty(zipfile),
    zipfile = ...
        'http://www.cs.tut.fi/~gomezher/tutorial-circuits-data.zip';
end

if nsamples > MAX_SAMPLES,
    error('The number of data samples must be smaller or equal to %d', ...
        MAX_SAMPLES);
end

if ntrials > MAX_TRIALS,
    error('The number of data trials must be smaller or equal to %d', ...
        MAX_TRIALS);
end

files = unzip(zipfile);

id = {};
pset = {};
count = 0;
for i = 1:length(files),
    % first circuit
    if length(files{i})>11 && strcmpi(files{i}(1:11),'inputMGScan'),
        index0 = strfind(files{i},'_');
        index0 = index0(end);
        index1 = strfind(files{i},'.');
        index1 = index1(end);
        idstr = files{i}((index0+1):index1-1);
        id = [id;{idstr}];
        tmp = load(files{i});
        pset = [pset,{squeeze(tmp(1:nsamples))'}];
        count = count + 1;
    end
    if count >= ntrials,
        break;
    end
end

% Initialize the output pointset
pset = [pset;cell(1,ntrials)];
count = 0;
for i = 1:length(files),
    % second circuit
    if length(files{i})>12 && strcmpi(files{i}(1:12),'outputMGScan'),
        index0 = strfind(files{i},'_');
        index0 = index0(end);
        index1 = strfind(files{i},'.');
        index1 = index1(end);
        idstr = files{i}((index0+1):index1-1);
        idstr = idstr(1:end-1);
        tmp = load(files{i});
        if any(ismember(id,idstr)),
            pset(2,ismember(id,idstr)) = {squeeze(tmp(1:nsamples))'};
            count = count + 1;
        end
    end
    if count >= ntrials,
        break;
    end
end

for i = 1:length(files),
    delete(files{i});
end









