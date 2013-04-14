function recurrence_plot(pointSet, varargin)

import([tim_package, '.*']);

concept_check(nargin, 'inputs', 1);
concept_check(nargout, 'outputs', 0 : 1);

% Optional input arguments
deviation = 8;
eval(process_options({'deviation'}, varargin));

% Compute the squared distance matrix.
distanceMatrix2 = tim.distance_matrix(pointSet);
distanceMatrix2 = exp(-distanceMatrix2 / (2 * deviation^2));

% Draw the Gaussian recurrence plot.
iptsetpref('ImshowAxesVisible', 'on');
figure;
imshow(distanceMatrix2);
axis xy;
title('Recurrence plot');
xlabel('j');
ylabel('i');
