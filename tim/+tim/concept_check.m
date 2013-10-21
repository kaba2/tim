% CONCEPT_CHECK
% Checks that a variable conforms to a given concept.
%
% concept_check(X, concept, parameter)
%
% where
%
% X is some variable.
%
% CONCEPT is the name of the concept. The available
% concepts are listed below. 
%
% PARAMETER is a parameter required by some, but
% not all concepts. In the latter case it can be
% omitted.
%
% Concepts
% --------
% 
% INPUTS ('inputs')
%
% The number of output arguments of a function.
%
% Parameter: An arbitrary-dimensional integer array
% denoting the set of legal number of input parameters.
%
% OUTPUTS ('outputs')
%
% The number of output arguments of a function.
%
% Parameter: An arbitrary-dimensional integer array
% denoting the set of legal number of output parameters.
%
% FILTER ('filter')
%
% A weighting function for temporal estimators. 
% An arbitrary-dimensional real array, whose linearization
% contains the weighting coefficients. Must have an odd
% number of elements.
%
% SIGNALSET ('tim.signal_set')
%
% An arbitrary-dimensional cell-array of signals.
% The dimensions of the trials of a given signal
% must match. The dimensions between signals may 
% or may not match. No signal must have dimension
% greater than 32.
%
% K ('k')
%
% The number of nearest neighbors to use in estimators.
%
% TIMEWINDOWRADIUS ('timeWindowRadius')
%
% The radius of the time window.
% An integer.

% Description: Checks that a variable conforms to a given concept.
% Documentation: tim_matlab_impl.txt

function concept_check(X, concept, parameter)

handled = 0;

maxDimension = 32;

if strcmp(concept, 'inputs')
    if isempty(find(parameter(:) == X, 1))
        if X < parameter
            error('Not enough input arguments.');
        end
        if X > parameter
            if mod(X - parameter, 2) ~= 0
                error(['The optional input arguments must be ', ...
                    'given in key-value pairs.']);
            end
        end
    end
    
    handled = 1;
end

if strcmp(concept, 'outputs')
    if isempty(find(parameter(:) == X, 1))
        if X < min(parameter(:))
            error('Not enough output arguments.');
        end
        if X > max(parameter(:))
            error('Too many output arguments.');
        end
        error(['The number of output arguments must be one of ', ...
            num2str(parameter(:)'), '.']);
    end
    handled = 1;
end

if ~handled
    error(['There is no concept check named "', concept, '".']);
end
