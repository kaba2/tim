% COMPUTE_DIRECTION_SET
% Computes the normalized differences between sub-sequent points.
%
% compute_direction_set(pointSet, 'key', value, ...)
%
% POINTSET is a real (d x n)-array where each column contains
% a d-dimensional point.

function directionSet = compute_direction_set(pointSet, varargin)

% Optional input arguments.
eval(tim.process_options(...
    {}, ...
    varargin));

[d, n] = size(pointSet);

directionSet = zeros(d, n);
directionSet(:, 1 : n - 1) = diff(pointSet, 1, 2);
directionSet(:, n) = directionSet(:, n - 1);

normSet = sqrt(sum(directionSet.^2));
directionSet = directionSet ./ (ones(d, 1) * normSet);
