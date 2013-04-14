% LORENZ
% Generates the Lorenz point-set.
%
% pointSet = lorenz();
% pointSet = lorenz('key', value, ...);
%
% Optional arguments
% ------------------
%
% SIGMA ('sigma') is a real number; a constant of the Lorenz system.
% Default: 10
%
% R ('r') is a real number; a constant of the Lorenz system.
% Default: 28
%
% B ('b') is a real number; a constant of the Lorenz system.
% Default: 8 / 3
%
% TMAX ('tMax') is a real number which contains the maximum
% time to cover. The time range will be [0, tMax].
% Default: 40
%
% N ('n') is an integer specifying the number of points to return.
% The points will be uniformly distributed in time.
% Default: ceil(100 * tMax)
%
% P0 ('p0') is the initial point representing time zero.
% Default: [1, 0, 0]
%
% The Lorenz point-set is the solution of the Lorenz
% ordinary differential equation
%
%     x' = sigma(y - x)
%     y' = rx - y - xz
%     z' = -bz + xy.

% Description: Generates the Lorenz point-set

function pointSet = lorenz(varargin)

import([tim_package, '.*']);

concept_check(nargin, 'inputs', 0);
concept_check(nargout, 'outputs', 0 : 1);

% Optional input arguments.
sigma = 10;
r = 28;
b = 8 / 3;
p0 = [1, 0, 0];
tMax = 40;
n = ceil(100 * tMax);
eval(process_options(...
    {'sigma', 'r', 'b', 'p0', 'tMax', 'n'}, ...
    varargin));

f = @(t, p) [...
    sigma * (p(2) - p(1)); ...
    r * p(1) - p(2) - p(1) * p(3); ...
    -b * p(3) + p(1) * p(2)];

tSet = linspace(0, tMax, n);
[ignore, pointSet] = ode45(f, tSet, p0);

pointSet = pointSet';
