% CHECK_SIGNALSET
% Checks a signal-set for basic validity.
%
% check_signalset(X)
%
% where
%
% X is an 2-dimensional array of signals.

% Description: Checks a signal-set for basic validity
% Documentation: tim_matlab_matlab.txt

function check_signalset(X)

if ~iscell(X)
    error('A signal-set is not a cell-array.');
end

% The limit for the dimension is arbitrary, but
% protects for the case when the user accidentally
% passes the transpose of the intended data.

maxDimension = 32;

signals = numel(X);

if signals == 0
	error('Signal-set is empty.');
end

dimension = size(X{1}, 1);

for i = 1 : signals
    if ~isa(X{i}, 'double')
        error('Some signal is not of type double.');
    end

    if size(X{i}, 1) ~= dimension
        error(['The dimensions of the trials do not match.']);
    end

    if size(X{i}, 1) > maxDimension
        error(['Some signal has dimension greater than ', ...
            int2str(maxDimension), '.']);
    end
end

end
