% POINCARE_SECTION
% Poincare section of a time-series.
%
% sectionSet = poincare_section(pointSet, 'key', value, ...)
%
% where
%
% POINTSET is a real (d x n)-array; a signal.
%
% Returned values
% ---------------
%
% SECTIONSET is a real ((d-1) x m)-array, which contains the interpolated
% crossings of the poincare plane.
%
% Optional arguments
% ------------------
%
% NORMAL ('normal') is a d-dimensional normal vector of the poincare plane.
% Default: [1; 0; 0]
%
% POSITION ('position') is a d-dimensional point through which the
% plane should pass.
% Default: [0; 0; 0]
%
% DIRECTION ('direction') is an integer in {-1, 0, 1}, which denotes
% the direction in which the poincare section is taken.
%  0: All intersections.
% +1: Only those intersections which agree with the normal.
% -1: Only those intersections which do not agree with the normal.
%
% More
% ----
%
% Let f : R^d --> {-1, 0, 1} be given by
% 
%     f(x) = sgn(dot(normal, (x - position))).
%
% A _crossing index_ of the poincare plane is an index i such that
%
%     i is a crossing <=> f(x(i)) f(x(i + 1)) < 0.
%
% Let 
%
%     p_i : R --> R^D: p_i(t) = x(i) + t (x(i + 1) - x(i)).
%
% Given a crossing index i, a _crossing point_ is a point
% p_i(t') in R^d defined by
%
%     dot(normal, (p_i(t') - position)) = 0
%     <=>
%     t' = dot(normal, position - x(i)) / dot(normal, x(i + 1) - x(i))
%
% i.e. the crossing point is interpolated between the points
% x(i) and x(i + 1) to yield a point on the poincare plane.

% Description: Poincare section of a time-series.

function sectionSet = poincare_section(pointSet, varargin)

import([tim_package, '.*']);

concept_check(nargin, 'inputs', 1);
concept_check(nargout, 'outputs', 0 : 1);

% Optional input arguments.
normal = [1; 0; 0];
position = [0; 0; 0];
direction = 0;
eval(process_options(...
    {'normal', 'position', 'direction'}, ...
    varargin));

% Ensure these are column-vectors.
normal = normal(:);
position = position(:);

d = size(pointSet, 1);
n = size(pointSet, 2);

% Generate an orthonormal basis on the plane
normal = normal / norm(normal);
planeBasis = null(normal');

sectionSet = zeros(d - 1, 0);
for i = 1 : n - 1
    currentSide = dot(normal, pointSet(:, i) - position);
    nextSide = dot(normal, pointSet(:, i + 1) - position);
    crossing = currentSide * nextSide < 0;
    rightDirection = (direction == 0 || ...
        (direction < 0) == (nextSide < 0));
    if crossing && rightDirection       
        t = dot(normal, position - pointSet(:, i)) / ...
            dot(normal, pointSet(:, i + 1) - pointSet(:, i));
        
        % Compute the crossing point by interpolation.
        crossingPoint = pointSet(:, i) + ...
            t * (pointSet(:, i + 1) - pointSet(:, i));
        
        % Project the point on the plane coordinates
        % (which are (d-1)-dimensional).
        sectionSet(:, end + 1) = ...
            planeBasis' * crossingPoint;
    end
end
