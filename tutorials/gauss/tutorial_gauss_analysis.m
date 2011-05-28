% TUTORIAL_GAUSS_ANALYSIS
% Estimates partial transfer entropy (PTE) for the toy experiment with three
% non-linearly coupled Gaussian processes
%
% results = tutorial_gauss_analysis(data)
%
% where
%
% DATA is a cell array generated using function gaussian_process.
% Alternatively, data can also be an string with the name of a .MAT file
% where the data is stored (in a variable named "data").
%
% AVG_FLAG is used to determine whether to use the "average estimators"
% (AVG_FLAG = 1) or the "ensemble estimators". Default 0.
%
% RESULTS is a struct containing the resuls of the analysis.
%
% See also: tutorial_gauss_figures, tutorial_gauss_data

% Description: Time-varying partial transfer entropy
% Documentation: tutorial_gauss.txt
% Author: German Gomez-Herrero

function results = tutorial_gauss_analysis(data)

IN_FILE = 'gaussian_process.mat';     % default data file
OUT_FILE = 'tutorial_gauss_analysis'; % results will be saved here
DIM = 1;                              % embedding dimension
STEP = 1;
K = 20;              % number of nearest neighbors
WINDOW_RADIUS = 5;
LAG21 = 9;
LAG32 = 14;
LAG31 = 24;
FILTER_ORDER = 20;   % order of the post-processing moving average filter
ALPHA = 0.1;        % significance level

if nargin < 1 || isempty(data),
    s = load(IN_FILE);
    data = s.data;
elseif ischar(data),
    s = load(data);
    data = s.data;
end

% Delay-embed the three time-series
S1 = delay_embed(data(1,:), DIM, STEP);
S2 = delay_embed(data(2,:), DIM, STEP);
S3 = delay_embed(data(3,:), DIM, STEP);
clear data;

% Function handle to an information transfer estimator
estimator = @(x, lag) transfer_entropy_pt(x(1,:), x(2,:), x(3,:), x(4,:), ...
    WINDOW_RADIUS, 'yLag', lag, 'k', K);
estimator_smooth = @(x, lag) filter((1/FILTER_ORDER)*ones(1,FILTER_ORDER),1, ...
    estimator(x, lag),[],2);

% Information transfer in the direction 2<-1
W = delay_embed_future(S2, STEP);
pte21 = estimator_smooth([S2;S1;S3;W],LAG21);

% Significance threshold for flow 2<-1
pte21sig = permutation_test([S2;S1;S3;W], [1 2 1 1], ...
    @(x) estimator_smooth(x, LAG21), ALPHA);

%  1<-2
W = delay_embed_future(S1, STEP);
pte12 = estimator_smooth([S1;S2;S3;W],0);
pte12sig = permutation_test([S2;S1;S3;W], [1 2 1 1], ...
    @(x) estimator_smooth(x, 0), ALPHA);

% 3<-2
W = delay_embed_future(S3, STEP);
pte32 = estimator_smooth([S3;S2;S1;W],LAG32);
pte32sig = permutation_test([S3;S2;S1;W], [1 2 1 1], ...
    @(x) estimator_smooth(x, LAG32), ALPHA);

% 2<-3
W = delay_embed_future(S2, STEP);
pte23 = estimator_smooth([S2;S3;S1;W],0);
pte23sig = permutation_test([S2;S3;S1;W], [1 2 1 1], ...
    @(x) estimator_smooth(x, 0), ALPHA);

% 1<-3
W = delay_embed_future(S1, STEP);
pte13 = estimator_smooth([S1;S3;S2;W],0);
pte13sig = permutation_test([S1;S3;S2;W], [1 2 1 1], ...
    @(x) estimator_smooth(x, 0), ALPHA);

% 3<-1
W = delay_embed_future(S3, STEP);
pte31 = estimator_smooth([S3;S1;S2;W],LAG31);
pte31sig = permutation_test([S3;S1;S2;W], [1 2 1 1], ...
    @(x) estimator_smooth(x, LAG31), ALPHA);

results = struct('pte12', pte12, 'pte12sig', pte12sig, ...
    'pte21', pte21, 'pte21sig', pte21sig, ...
    'pte23', pte23, 'pte23sig', pte23sig, ...
    'pte32', pte32, 'pte32sig', pte32sig, ...
    'pte13', pte13, 'pte13sig', pte13sig, ...
    'pte31', pte31, 'pte31sig', pte31sig, ...
    'windowRadius', WINDOW_RADIUS, 'filterOrder', FILTER_ORDER, ...
    'k', K, 'alpha', ALPHA, 'step', STEP, 'dim', DIM);


save([OUT_FILE '.mat'], 'results');
