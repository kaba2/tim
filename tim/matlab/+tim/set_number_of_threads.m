% SET_NUMBER_OF_THREADS
% Sets the number of threads to use for computation.
%
% set_number_of_threads(threads)
%
% where
%
% THREADS is a positive integer denoting the number of
% threads to use for computation.
%
% Type 'help tim' for more documentation.

% Description: Sets the number of threads to use for computation
% Documentation: support_functions.txt

function set_number_of_threads(threads)

import([tim_package, '.*']);

concept_check(nargin, 'inputs', 1);

tim_matlab('set_number_of_threads', threads);
