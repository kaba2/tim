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
% TSET ('tSet') is a real-array, whose linearization contains 
% the time-points on which to evaluate the Lorenz point-set.
% Default: 0 : 0.01 : 100
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
tSet = 0 : 0.01 : 100;
eval(process_options(...
    {'sigma', 'r', 'b', 'p0', 'tSet'}, ...
    varargin));

f = @(t, p) [...
    sigma * (p(2) - p(1)); ...
    r * p(1) - p(2) - p(1) * p(3); ...
    -b * p(3) + p(1) * p(2)];

[ignore, pointSet] = ode45(f, tSet, p0);

pointSet = pointSet';
