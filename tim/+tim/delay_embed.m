% DELAY_EMBED 
% Delay-embeds a signal.
%
% Y = delay_embed(X, k, tDelta)
%
% where
%
% K is the "embedding factor".
%
% TDELTA is the embedding delay.

% Description: Delay-embedding
% Documentation: delay_embed.txt

function Y = delay_embed(X, k, tDelta)

if nargin < 3 || isempty(tDelta), tDelta = 1; end
if nargin < 2 || isempty(k) || isempty(X),
    error('Not enough input arguments');
end

% Deal with the case of multiple input signals.

if iscell(X),
    Y = cell(size(X));
    for i = 1:length(X)
        Y{i} = delay_embed(X{i}, k, tDelta);        
    end
    return;
end

if k < 0,
    error('k must be a non-negative integer.');
end

if tDelta < 1,
    error('tDelta must be a positive integer.');
end

% use this line when TIM estimators admit NaNs
X = [nan(size(X, 1), (k - 1) * tDelta), X];
% for the moment we use circular shifting
%tmp = circshift(X, [0 (k-1)*tDelta]);
%X = [tmp(:,1:(k-1)*tDelta) X];

n = size(X, 1);

embedDimension = k * n;
embedSampleWidth = (k - 1) * tDelta + 1;
embedSamples = size(X, 2) - embedSampleWidth + 1;

Y = zeros(embedDimension, embedSamples);

for j = 1 : k
   s = tDelta * (j - 1) + 1; 
   Y((j - 1) * n + 1 : j * n, :) = X(:, s : (s + embedSamples - 1));
end
