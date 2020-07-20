% COUPLED_LORENZ_SYSTEMS
% Generates two coupled Lorenz point-sets.
%
% [aSet, bSet] = coupled_lorenz_systems();
% [aSet, bSet] = coupled_lorenz_systems('key', value, ...);
%
% Return values
% -------------
%
% ASET and BSET are real (3 x n)-matrices each containing n 
% 3-dimensional points as columns.
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
% Default: [0, 1, 0, 0, 1, 0]
%
% ALPHA ('alpha') is a real number which specifies how strongly
% the second Lorenz system drives the first one.
% Default: 0
%
% BETA ('beta') is a real number which specifies how strongly
% the first Lorenz system drives the second one.
% Default: 0.5
%
% Additional information
% ----------------------
%
% The coupled Lorenz point-sets are the solution to the following
% ordinary differential equation:
%
%     x_1' = sigma y_1 - sigma(x_1 - alpha x_2)
%     y_1' = r x_1 - y_1 - x_1 z_1
%     z_1' = -b z_1 + x_1 y_1.
%
%     x_2' = sigma y_2 - sigma(x_2 - beta x_1)
%     y_2' = r x_2 - y_2 - x_2 z_2
%     z_2' = -b z_2 + x_2 y_2.

% Description: Generates two coupled Lorenz point-sets

function [aSet, bSet] = coupled_lorenz_systems(varargin)

import([tim_package, '.*']);

concept_check(nargin, 'inputs', 0);
concept_check(nargout, 'outputs', 0 : 2);

% Optional input arguments.
sigma = 10;
r = 28;
b = 8 / 3;
p0 = [0, 1, 0, 0, 1, 0];
tMax = 40;
n = ceil(100 * tMax);
alpha = 0;
beta = 0.5;
eval(process_options(...
    {'sigma', 'r', 'b', 'p0', 'tMax', 'n', ...
	'alpha', 'beta'}, ...
    varargin));

pastelmatlab.concept_check(...
	sigma, 'real', ...
	r, 'real',  ...
	b, 'real', ...
	tMax, 'real', ...
	n, 'integer', ...
	alpha, 'real', ...
	beta, 'real')

f = @(t, p) [...
    sigma * p(2) - sigma * (p(1) - alpha * p(4)); ...
    r * p(1) - p(2) - p(1) * p(3); ...
    -b * p(3) + p(1) * p(2); ...
    sigma * p(5) - sigma * (p(4) - alpha * p(1)); ...
    r * p(4) - p(5) - p(4) * p(6); ...
    -b * p(6) + p(4) * p(5)];

tSet = linspace(0, tMax, n);
[ignore, pointSet] = ode45(f, tSet, p0);

pointSet = pointSet';

aSet = pointSet(1 : 3, :);
bSet = pointSet(4 : 6, :);

