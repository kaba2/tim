% AUTO_MI
% Computes mutual information with lagged versions of itself.
%
% miSet = auto_mi(pointSet)
% miSet = auto_mi(pointSet, 'key', value, ...)
%
% where
%
% POINTSET is a real (d x n)-matrix whose columns contain
% d-dimensional points.
%
% MISET is a real (1 x n)-matrix where MISET(i) is the mutual-information
% between POINTSET and POINTSET lagged with LAGSET(i).
%
% Optional arguments
% ------------------
%
% LAGSET ('lagSet') is an integer array whose linearization 
% contains p lags. Default: [0 : 63]
%
% Type 'help tim' for more documentation.

function miSet = auto_mi(pointSet, varargin)

import([tim_package, '.*']);

concept_check(nargin, 'inputs', 1);
concept_check(nargout, 'outputs', 0 : 1);

% Optional input arguments.
lagSet = 0 : 63;
eval(tim.process_options(...
    {'lagSet'}, ...
    varargin));

% Compute the auto mutual information for the
% given lags.
miSet = tim.mutual_information(...
    pointSet, pointSet, ...
    'yLag', lagSet)';

