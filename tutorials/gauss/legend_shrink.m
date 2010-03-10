% LEGEND_SHRINK 
% Shrinks the size of the legend
%
% legend_shrink(legend_handle)
%
% where LEGEND_HANDLE is a handle to a legend object. 

% Description: Shrink the size of the legend of a figure
% Documentation: tutorial_gauss.txt
% Author: German Gomez-Herrero

function y = legend_shrink(h, factor, factor2, factor3)

if nargin < 4 || isempty(factor3), factor3 = .95; end
if nargin < 3 || isempty(factor2), factor2 = .8; end
if nargin < 2 || isempty(factor), factor = .6; end

children = get(h, 'Children');

for i = 1:length(children)
   type = get(children(i), 'type');    
   switch lower(type),
       case 'line'
           % modify line lengths
           xdata = get(children(i), 'XData');           
           if length(xdata) > 1,
               xdata(2) = xdata(2)*factor;
           end
           set(children(i), 'XData', xdata);
           
           
       case 'hggroup',
           % area line+tag
           morechildren = get(children(i), 'Children');
           for j = 1:length(morechildren),
               if strcmpi(get(morechildren(j), 'type'), 'patch')
                   vertices = get(morechildren(j), 'Vertices');  
                   vertices(2,2) = vertices(1,2)+(vertices(2,2)-vertices(1,2))*factor2;
                   vertices(3,2) = vertices(4,2)+(vertices(3,2)-vertices(4,2))*factor2;
                   vertices(3,1) = vertices(3,1)*factor;
                   vertices(4,1) = vertices(4,1)*factor;
                   set(morechildren(j), 'Vertices', vertices);
               end
           end
           
           
       case 'text',
           position = get(children(i), 'Position');
           position(1) = position(1)*factor;
           set(children(i), 'Position', position);           
           
           
       otherwise,
           
   end
end
position = get(h, 'Position');
%position(3) = position(3) * factor;
%outerPosition = get(h, 'OuterPosition');
%outerPosition(3) = outerPosition(3) * factor;
dataAspectRatio = get(h, 'DataAspectRatio');
dataAspectRatio(2) = position(3)/position(4) * factor3;
set(h, 'DataAspectRatio', dataAspectRatio);
%set(h, 'position', position);