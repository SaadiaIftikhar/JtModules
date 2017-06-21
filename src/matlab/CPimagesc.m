function h = CPimagesc(Image,handles,matCLim)

% CellProfiler is distributed under the GNU General Public License.
% See the accompanying file LICENSE for details.
%
% Developed by the Whitehead Institute for Biomedical Research.
% Copyright 2003,2004,2005.
%
% Authors:
%   Anne E. Carpenter
%   Thouis Ray Jones
%   In Han Kang
%   Ola Friman
%   Steve Lowe
%   Joo Han Chang
%   Colin Clarke
%   Mike Lamprecht
%   Peter Swire
%   Rodrigo Ipince
%   Vicky Lay
%   Jun Liu
%   Chris Gang
%
% Website: http://www.cellprofiler.org
%
% $Revision: 2802 $

%%% Displays the image.
% [BS-HACK, made intensity-range quantile based to actually show something!]
% h = imagesc(Image,[quantile(Image(:),0.02), quantile(Image(:),0.98)]);
% [NB-HACK], add posibillity of passing the limits for the imagesc function
if islogical(Image)
    h = imagesc(Image);
else
    if nargin<3
        matCLim = [quantile(Image(:),0.001) quantile(Image(:),0.999)];
    end
    
    if matCLim(1)~=matCLim(2)
        h = imagesc(Image,[matCLim(1), matCLim(2)]);
    else
        h = imagesc(Image);
    end
end


%%% Embeds the Image tool submenu so that it appears when the user clicks on the image.
set(h,'ButtonDownFcn','CPimagetool');

%%% Sets the user's preference for font size, which should affect tick
%%% labels and current and future titles.
set(gca,'fontsize',handles.Preferences.FontSize)

%%% Applies the user's choice for colormap.
if ndims(Image) == 2
    colormap(handles.Preferences.IntensityColorMap);
end