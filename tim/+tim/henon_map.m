% HENON
% Generates the Henon point-set.
%
% pointSet = henon();
% pointSet = henon('key', value, ...);
%
% Optional arguments
% ------------------
%
% A ('a') is a real number; a constant of the Henon map.
% Default: 1.4
%
% B ('b') is a real number; a constant of the Henon map.
% Default: 0.3
%
% P0 ('p0') is the initial point representing time zero.
% Default: [0, 0]
%
% N ('n') is the number of points to generate.
% Default: 1000
%
% The Henon map is an iterated function given by
%
%     x_{i+1} = a - x_i^2 + b y_i,
%     y_{i+1} = x_i,
%
% where
%
%     x_0 = p0(1),
%     y_0 = p0(2).

% Description: Generates the Lorenz point-set

function pointSet = henon_map(varargin)

import([tim_package, '.*']);

concept_check(nargin, 'inputs', 0);
concept_check(nargout, 'outputs', 0 : 1);

% Optional input arguments.
a = 1.4;
b = 0.3;
p0 = [0, 0];
n = 1000;
eval(process_options(...
    {'a', 'b', 'p0', 'n'}, ...
    varargin));

pointSet = zeros(2, n);
x = p0(1);
y = p0(2);
for i = 1 : n
    pointSet(:, i) = [x; y];
    yPrev = y;
    y = x;
    x = a - x^2 + b * yPrev;
end
