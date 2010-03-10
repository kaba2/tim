% CHECK
% Checks a variable for validity.
%
% check(X, check)
%
% where
%
% X is some variable.
%
% CHECK is one of:
% 'filter'
% 'signalSet'
% 'k'
% 'timeWindowRadius'
% 'threads'

% Description: Checks a filter for basic validity
% Documentation: tim_matlab_matlab.txt

function check_variable(X, check)

handled = 0;

if strcmp(check, 'filter')
	if ~isa(X, 'double')
		error('FILTER must be a real array');
	end

	if mod(numel(X), 2) == 0
		error('FILTER must have odd number of elements.');
	end	

	if sum(X(:)) == 0
		error('FILTER must not sum to 0');
    end	
    handled = 1;
end

if strcmp(check, 'signalSet')
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
    
    handled = 1;
end

if strcmp(check, 'k')
	if size(X, 1) ~= 1 || ...
	   size(X, 2) ~= 1
		error('K must be a scalar integer.');
	end

	if X < 1
		error('K must be at least 1.');
    end
    
    handled = 1;
end

if strcmp(check, 'timeWindowRadius')
	if size(X, 1) ~= 1 || ...
	   size(X, 2) ~= 1
		error('TIMEWINDOWRADIUS must be a scalar integer.');
	end

	if timeWindowRadius < 0
		error('TIMEWINDOWRADIUS must be non-negative');
    end
    handled = 1;
end

if strcmp(check, 'threads')
    if size(threads, 1) ~= 1 || ...
       size(threads, 2) ~= 1
        error('THREADS must be a scalar integer.');
    end

    if threads < 1
        error('THREADS must be at least 1.');
    end
    handled = 1;
end

if ~handled
    error(['No such check ', check, '.']);
end