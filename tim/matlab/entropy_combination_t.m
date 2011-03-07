% ENTROPY_COMBINATION_T
% A temporal entropy combination estimate from samples.
%
% I = entropy_combination_t(signalSet, rangeSet, 
%       timeWindowRadius, lagSet, k, filter)
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
% LAGSET is an arbitrary-dimensional cell-array whose linearization 
% contains p arrays of lags to apply to each signal. Each array of lags
% is either a scalar, or has L elements, where L is the maximum number of 
% elements among the arrays in LAGSET. An array of lags is handled by its
% linearization. If an array of lags is a scalar, it is extended to an
% array with L elements with the scalar as its elements.
% Default: a (p x 1) cell-array of scalar zeros.
%
% Type 'help tim' for more documentation.

% Description: Temporal entropy combination estimation
% Documentation: entropy_combination.txt

function I = entropy_combination_t(signalSet, rangeSet, ...
    timeWindowRadius, lagSet, k, filter)

check(nargin, 'inputs', 3 : 6);
check(nargout, 'outputs', 0 : 1);

if nargin < 4
	lagSet = num2cell(zeros(size(signalSet, 1), 1));
end

if nargin < 5
    k = 1;
end

if nargin < 6
	filter = [1];
end

signals = size(signalSet, 1);

for i = 1 : signals
    check(signalSet(i, :), 'signalSet');
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

check(k, 'k')
check(filter, 'filter');

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

