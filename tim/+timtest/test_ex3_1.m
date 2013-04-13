% Non-linear time series analysis
% Exercise 3.1

clear all;
close all;

pointSet = tim.henon('n', 2000);

% The first minimum of auto mutual information fails to
% predict a good embedding lag for this data set.
% The reason is that the Henon map lacks continuity, an 
% assumption of the automated method. The best embedding lag 
% is empirically found to be 1 (just try different choices).
% I suspect this is a general property of an iterated
% function; every step is important.

% The auto mutual information goes to zero almost immediately; 
% again I suspect this is a general property of an iterated function.

timtest.visualize_embedding(pointSet, 'tDelta', 1);

