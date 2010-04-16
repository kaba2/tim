% DELAY_EMBED_FUTURE
% Future of a delay-embedded signal.
%
% Y = delay_embed_future(X, k, dt)
%
% where
%
% X is a signal before delay-embedding.
%
% K is the "embedding factor".
%
% DT is the embedding delay.

% Description: Future of a delay-embedded signal
% Documentation: tim_matlab_matlab.txt

function Y = delay_embed_future(X, k, dt)

if nargin < 3 || isempty(dt), dt = 1; end
if nargin < 2 || isempty(k) || isempty(X),
    k = 1;
end

% Deal with the case of multiple input signals.

if iscell(X),
    Y = cell(size(X));
    for i = 1:length(X)
        Y{i} = delay_embed_future(X{i}, k, dt);        
    end
    return;
end

if k < 0,
    error('k must be a non-negative integer.');
end

if dt < 1,
    error('dt must be a positive integer.');
end

futureShift = dt * k;

Y = X(:, futureShift + 1 : end);
