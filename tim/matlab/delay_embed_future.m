% DELAY_EMBED_FUTURE
% Future of a delay-embedded signal.
%
% Y = delay_embed_future(X, k, t0, dt)
%
% where
%
% X is a signal before delay-embedding.
%
% K is the "embedding factor".
%
% T0 is the embedding t0.
%
% DT is the embedding delay.

% Description: Future of a delay-embedded signal
% Documentation: tim_matlab_matlab.txt

function Y = delay_embed_future(X, k, t0, dt)

if nargin < 4 || isempty(dt), dt = 1; end
if nargin < 3 || isempty(t0), t0 = 0; end
if nargin < 2 || isempty(k) || isempty(X),
    error('Not enough input arguments');
end

% Deal with the case of multiple input signals.

if iscell(X),
    Y = cell(size(X));
    for i = 1:length(X)
        Y{i} = delay_embed_future(X{i}, k, t0, dt);        
    end
    return;
end

if k < 0,
    error('k must be a non-negative integer.');
end

if t0 < 0,
    error('t0 must be a non-negative integer.');
end

if dt < 1,
    error('dt must be a positive integer.');
end

future_shift = t0 + dt * k;
Y = delay_embed(X, 1, future_shift);