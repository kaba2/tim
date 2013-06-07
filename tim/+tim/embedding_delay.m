% EMBEDDING_DELAY
% Finds a good embedding delay
%
% dt = embedding_delay(pointSet, 'key', value, ...)
%
% where
%
% POINTSET is a real (d x n)-matrix which contains n d-dimensional
% points.
%
% Optional arguments
% ------------------
%
% MAXLAG ('maxLag') is a positive integer which specifies the 
% maximum lag to inspect.
% Default: size(pointSet, 2)
%
% Additional information
% ----------------------
%
% A good embedding delay is given by the first minimum of the
% auto mutual information. If the signal is sampled densely enough, 
% temporally close samples are related simply because of continuity. 
% Increasing the lag decreases the continuity-caused dependency
% between delayed samples. On the other hand, at some point the
% dependency starts to rise again, since a chaotic system usually
% has cyclic-like properties. Overall, the dependency will lower,
% since temporally distant samples have less and less dependency to
% each other, because of the chaotic pseudo-randomness. This is
% why one should choose the first minimum, rather than subsequent
% minima.

function dt = embedding_delay(pointSet, varargin)

import([tim_package, '.*']);

concept_check(nargin, 'inputs', 1);
concept_check(nargout, 'outputs', 0 : 1);

[d, n] = size(pointSet);

% Optional input arguments.
maxLag = n;
eval(process_options(...
    {'maxLag'}, ...
    varargin));

prevMi = Inf;
yLag = 0;
while yLag <= maxLag;
	mi = tim.mutual_information(...
	    pointSet, pointSet, ...
	    'yLag', yLag + 1);
	
	if prevMi < mi
		% We are at a local minimum.
		% Return the current yLag  
		% as the embedding delay.
		break;
	end

	prevMi = mi;
	yLag = yLag + 1;
end

dt = yLag;
