% CHECK
% Checks a variable for validity.
%
% check(X, checkName)
%
% where
%
% X is some variable.
%
% CHECK is one of:
% 'inputs'
% 'outputs'
% 'filter'
% 'signalSet'
% 'k'
% 'timeWindowRadius'

% Description: Checks a variable for validity
% Documentation: tim_matlab_impl.txt

function check(X, checkName, checkParam)

handled = 0;

maxDimension = 32;

if strcmp(checkName, 'inputs')
    if isempty(find(checkParam(:) == X, 1))
        if X < min(checkParam(:))
            error('Not enough input arguments.');
        end
        if X > max(checkParam(:))
            error('Too many input arguments.');
        end
        error(['The number of input arguments must be one of ', ...
            num2str(checkParam(:)'), '.']);
    end
    handled = 1;
end

if strcmp(checkName, 'outputs')
    if isempty(find(checkParam(:) == X, 1))
        if X < min(checkParam(:))
            error('Not enough output arguments.');
        end
        if X > max(checkParam(:))
            error('Too many output arguments.');
        end
        error(['The number of output arguments must be one of ', ...
            num2str(checkParam(:)'), '.']);
    end
    handled = 1;
end

if strcmp(checkName, 'filter')
	if ~isa(X, 'double')
		error('FILTER must be a real array');
	end

	if mod(numel(X), 2) == 0
		error('FILTER must have an odd number of elements.');
	end	

	if sum(X(:)) == 0
		error('FILTER must not sum to 0');
    end	
    handled = 1;
end

if strcmp(checkName, 'signalSet')
	if ~iscell(X)
		error('A signal-set is not a cell-array.');
	end

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

if strcmp(checkName, 'k')
	if size(X, 1) ~= 1 || ...
	   size(X, 2) ~= 1
		error('K must be a scalar integer.');
	end

	if X < 1
		error('K must be at least 1.');
    end
    
    handled = 1;
end

if strcmp(checkName, 'timeWindowRadius')
	if size(X, 1) ~= 1 || ...
	   size(X, 2) ~= 1
		error('TIMEWINDOWRADIUS must be a scalar integer.');
	end

	if X < 0
		error('TIMEWINDOWRADIUS must be non-negative');
    end
    handled = 1;
end

if ~handled
    error(['There is no check named "', checkName, '".']);
end