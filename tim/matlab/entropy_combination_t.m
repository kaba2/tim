% ENTROPY_COMBINATION_T
% A temporal entropy combination estimate from samples.
%
% I = entropy_combination_t(signalSet, rangeSet)
% I = entropy_combination_t(signalSet, rangeSet, 'key', value, ...)
%
% where 
%
% SIGNALSET is a 2-dimensional (p x q) cell-array 
% containing q trials of p signals.
%
% RANGESET is a 2-dimensional (m x 3) array which on each row 
% contains an integer triple of the form (a, b, s). Each triple
% denotes the rows (a : b) in SIGNALSET. Those signals are concatenated
% to form a higher-dimensional signal for which differential entropy
% is (implicitly) computed in the entropy combination. The s denotes
% the factor by which the differential entropy is multiplied before
% summing to the end-result.
%
% Optional input arguments in 'key'-value pairs:
%
% TIMEWINDOWRADIUS ('timeWindowRadius') is an integer that determines the
% temporal radius (in samples) around each point that will be used by the
% estimator. Default: ceil(200/q)
%
% LAGSET ('lagSet') is an arbitrary-dimensional cell-array whose 
% linearization contains p arrays of lags to apply to each signal. Each 
% array of lags is either a scalar, or has L elements, where L is the 
% maximum number of elements among the arrays in LAGSET. An array of lags
% is handled by its linearization. If an array of lags is a scalar, it is
% extended to an array with L elements with the scalar as its elements.
% Default: a (p x 1) cell-array of scalar zeros.
%
% K ('k') is an integer which denotes the number of nearest neighbors to be
% used by the estimator.
%
% FILTER ('filter') is an arbitrary-dimensional real-array, whose
% linearization contains temporal weighting coefficients. 
% Default: 1 (i.e. no temporal weighting is performed)
%
% Type 'help tim' for more documentation.

% Description: Temporal entropy combination estimation
% Documentation: entropy_combination.txt

function I = entropy_combination_t(signalSet, rangeSet, varargin)

pkgname = regexpi(mfilename('fullpath'), ['+(TIM_.*)' filesep mfilename], ...
    'tokens', 'once');
import([pkgname{1} '.tim_matlab']);
import([pkgname{1} '.concept_check']);
import([pkgname{1} '.process_options']);
import([pkgname{1} '.compute_lagarray']);

concept_check(nargin, 'inputs', 2);
concept_check(nargout, 'outputs', 0 : 1);

keySet = {'timeWindowRadius', 'lagSet', 'k', 'filter'};

% Default values of the optional input arguments
q = size(signalSet, 2);
timeWindowRadius = ceil(200/q);
lagSet = num2cell(zeros(size(signalSet, 1), 1));
k = 1;
filter = 1;

eval(process_options(keySet, varargin));

signals = size(signalSet, 1);

for i = 1 : signals
    concept_check(signalSet(i, :), 'signalSet');
end

marginals = size(rangeSet, 1);

if marginals == 0
    error('RANGESET is empty.');
end

if size(rangeSet, 2) ~=3
    error('The width of RANGESET must be 3');
end

for i = 1 : marginals
    if rangeSet(i, 1) > rangeSet(i, 2)
        error('For each RANGESET triple (a, b, s) must hold a <= b.');
    end
    if rangeSet(i, 1) < 1
        error('There is a RANGESET triple (a, b, s) with a < 1.');
    end
    if rangeSet(i, 2) > signals
        error('There is a RANGESET triple (a, b, s) with b > signals.');
    end
end

if numel(lagSet) ~= signals
	error(['LAGSET must contain the same number of elements as ', ...
	'there are signals in SIGNALSET.']);
end

lagArray = compute_lagarray(lagSet);

concept_check(k, 'k')
concept_check(filter, 'filter');

lags = size(lagArray, 2);

estimateSet = cell(1, lags);

for i = 1 : lags
    estimateSet{i} = tim_matlab(...
        'entropy_combination_t', ...
        signalSet, rangeSet, timeWindowRadius, ...
        lagArray(:, i), k, filter(:));
end

maxSamples = 0;
for i = 1 : lags
    if numel(estimateSet{i}) > maxSamples
        maxSamples = numel(estimateSet{i});
    end
end

I = nan(lags, maxSamples);
for i = 1 : lags
    I(i, 1 : numel(estimateSet{i})) = estimateSet{i};
end

