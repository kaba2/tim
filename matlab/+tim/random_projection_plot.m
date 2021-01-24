% RANDOM_PROJECTION_PLOT
% Visualizes a time-series by random projections.
%
% spatial_plot(pointSet, 'key', value, ...)
%
% POINTSET is a real (d x n)-array where each column contains
% a d-dimensional point.
%
% Optional arguments
% ------------------
%
% DIMENSION ('dimension'): An integer specifying the dimension of
% the projection subspaces. Must be 1, 2, or 3.
% Default: 2
%
% PLOTS ('plots'): A pair of integers, where plots(1) gives the
% number of rows and plots(2) gives the number of columns in an
% array of plots to produce. Each plot will be produced as a 
% subplot of the same figure.
% Default: [2, 2]
%
% SPATIAL_PLOT ('spatial_plot'): A cell-array of key-value pairs
% which specifies the optional arguments that are to be passed on
% to the spatial_plot() function (which draws the plots).
% Default: {}

% Description: Random-projection plot

function random_projection_plot(pointSet, varargin)

n = size(pointSet, 2);

% Optional input arguments.
dimension = 2;
plots = [2, 2];
spatial_plot = {};
eval(tim.process_options(...
    {'dimension', 'plots', 'spatial_plot'}, ...
    varargin));

for i = 1 : plots(1)
    for j = 1 : plots(2)
        k = sub2ind(plots, i, j);
        subplot(plots(1), plots(2), k);
        tim.spatial_plot(...
            pastelmatlab.random_projection(pointSet, ...
            'dimension', dimension), ...
            spatial_plot{:});
    end
end
