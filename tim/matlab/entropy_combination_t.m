% ENTROPY_COMBINATION_T
% A temporal entropy combination estimate from samples.
%
% I = entropy_combination_t(signalSet, rangeSet, 
%       timeWindowRadius, lagSet, k, threads)
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
% TIMEWINDOWRADIUS determines the radius of the time-window inside which
% samples are taken into consideration to the mutual information
% estimate at time instant t. The time-window at time instant t
% is given by [t - timeWindowRadius, t + timeWindowRadius]. 
% This allows the estimate to be adaptive to temporal changes in mutual 
% information. If no such changes should happen, better accuracy can be 
% achieved by either setting 'timeWindowRadius' to the number of samples 
% or using the entropy_combination() function instead.
%
% LAGSET is an arbitrary-dimensional array whose linearization contains
% p integers of lags to apply to each signal. Default: a (p x 1) array 
% of zeros.
%
% K determines which k:th nearest neighbor the algorithm
% uses for estimation. Default 1.
%
% THREADS determines the number of threads to use for parallelization.
% To fully take advantage of multiple cores in your machine, set this
% to the number of cores in your machine. Note however that this makes 
% your computer unresponsive to other tasks. When you need responsiveness, 
% spare one core for other work. Default maxNumCompThreads.
%
% Each signal is a real (m x n)-matrix that contains n samples of an
% m-dimensional signal. The dimensions of the signals need not coincide.
% However, the trials of a given signal must have equal dimension.
% If the number of samples varies with trials, the function uses 
% the minimum sample count among all the trials in SIGNALSET.

% Description: Temporal entropy combination estimation
% Documentation: tim_matlab_matlab.txt

function I = entropy_combination_t(signalSet, rangeSet, ...
    timeWindowRadius, lagSet, k, threads)

if nargin < 3
    error('Not enough input arguments.');
end

if nargin > 6
    error('Too many input arguments.');
end

if nargout > 1
    error('Too many output arguments.');
end

if nargin < 4
	lagSet = zeros(size(signalSet, 1), 1);
end

if nargin < 5
    k = 1;
end

if nargin < 6
    threads = maxNumCompThreads;
end

signals = size(signalSet, 1);

for i = 1 : signals
    check_signalset(signalSet(i, :));
end

if timeWindowRadius < 0
    error('TIMEWINDOWRADIUS must be non-negative');
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

lags = numel(lagSet);
if lags ~= signals
	error(['LAGSET must contain the same number of elements as ', ...
	'there are signals in SIGNALSET.']);
end

if size(k, 1) ~= 1 || ...
   size(k, 2) ~= 1
    error('K must be a scalar integer.');
end

if k < 1
    error('K must be at least 1.');
end

if size(threads, 1) ~= 1 || ...
   size(threads, 2) ~= 1
    error('THREADS must be a scalar integer.');
end

if threads < 1
    error('THREADS must be at least 1.');
end

I = tim_matlab('entropy_combination_t', ...
	signalSet, rangeSet, ...
    timeWindowRadius, lagSet, k, threads);

