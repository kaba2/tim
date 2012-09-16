% TIM 1.2.0
% Kalle Rutanen
% http://kaba.hilvi.org
% Copyright (c) 2009 - 2011
%
% This library is free software: you can redistribute it and/or modify
% it under the terms of the GNU Lesser General Public License as published 
% by the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This library is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
% GNU Lesser General Public License for more details.
% 
% You should have received a copy of the GNU Lesser General Public License
% along with this library. If not, see <http://www.gnu.org/licenses/>.

newEstimate = transfer_entropy_t(...
        signalSet(3, :), ...
        signalSet(1, :), ...
        signalSet(2, :), ...
        20, 'k', 10);

% estimateSet = tim_matlab(...
%         'entropy_combination_t', ...
%         signalSet, [1 2 1;2 3 1;2 2 -1], 20, ...
%         [0 0 0], 10, 1);
% 
