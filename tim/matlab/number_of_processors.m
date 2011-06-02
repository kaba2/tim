% NUMBER_OF_PROCESSORS
% Returns the number of processors in the machine
%
% processors = number_of_processors()
%
% where
%
% PROCESSORS is a positive integer. It is the number of 
% processors in the machine.
%
% Type 'help tim' for more documentation.

% Description: Number of processors in the machine
% Documentation: support_functions.txt

function processors = number_of_processors()

% Package initialization
eval(package_init(mfilename('fullpath')));

concept_check(nargin, 'inputs', 0);
concept_check(nargout, 'outputs', 1);

processors = tim_matlab('number_of_processors');
