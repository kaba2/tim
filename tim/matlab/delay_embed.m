% DELAY_EMBED 
% Delay-embeds a signal.
%
% Y = delay_embed(X, k, t0, dt)
%
% where
%
% K is the "embedding factor".
%
% T0 is the embedding shift.
%
% DT is the embedding delay.

% Description: Delay-embedding
% Documentation: tim_matlab_matlab.txt

function Y = delay_embed(X, k, t0, dt)

if nargin < 4 || isempty(dt), dt = 1; end
if nargin < 3 || isempty(t0), t0 = 0; end
if nargin < 2 || isempty(k) || isempty(X),
    error('Not enough input arguments');
end

% Deal with the case of multiple input signals.

if iscell(X),
    Y = cell(size(X));
    for i = 1:length(X)
        Y{i} = delay_embed(X{i}, k, t0, dt);        
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

n = size(X, 1);

embed_dimension = k * n;
embed_sample_width = (k - 1) * dt + 1;
extra_samples = t0 + embed_sample_width - 1;
X = [X X(:, 1 : extra_samples)];
embed_samples = (size(X, 2) - t0) - embed_sample_width + 1;

Y = zeros(embed_dimension, embed_samples);

for j = 1:k
   s = (t0 + dt * (j - 1) + 1); 
   Y((j - 1) * n + 1 : j * n, :) = X(:, s : (s + embed_samples - 1));
end
