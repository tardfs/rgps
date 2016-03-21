function nlib_skyplot(Az, El, sMap)

% This code is based on -*skyplot*- developed by Darius Plausinaitis 
% and Kristin Larson
% 
%--------------------------------------------------------------------------
%                           SoftGNSS v3.0
% 
% Copyright (C) Darius Plausinaitis and Kristin Larson
% Written by Darius Plausinaitis and Kristin Larson
%--------------------------------------------------------------------------
%This program is free software; you can redistribute it and/or
%modify it under the terms of the GNU General Public License
%as published by the Free Software Foundation; either version 2
%of the License, or (at your option) any later version.
%
%This program is distributed in the hope that it will be useful,
%but WITHOUT ANY WARRANTY; without even the implied warranty of
%MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%GNU General Public License for more details.
%
%You should have received a copy of the GNU General Public License
%along with this program; if not, write to the Free Software
%Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301,
%USA.
%--------------------------------------------------------------------------

%% Prepare axis ===========================================================
hAxis = gca ;

%--- Get x-axis text color so grid is in same color -----------------------
tc = get(hAxis, 'xcolor') ;

hold(hAxis, 'on') ;

%--- Plot white background ------------------------------------------------
rectangle('position', [-90, -90, 180, 180], ...
          'Curvature', [1 1], ...
          'facecolor', 'white', ...
          'edgecolor', tc);

%% Plot spokes ============================================================

%--- Find spoke angles ----------------------------------------------------
% Only 6 lines are needed to divide circle into 12 parts
th = (1:6) * 2*pi / 12;

%--- Convert spoke end point coordinate to Cartesian system ---------------
cst = cos(th); snt = sin(th);
cs = [cst; -cst];
sn = [snt; -snt];

%--- Plot the spoke lines -------------------------------------------------
line(90*sn, 90*cs, 'linestyle', ':', 'color', tc, 'linewidth', 0.5, ...
    'handlevisibility', 'off');

%% Annotate spokes in degrees =============================================
rt = 1.1 * 90;

for i = 1:max(size(th))

    %--- Write text in the first half of the plot -------------------------
    text(rt*snt(i), rt*cst(i), int2str(i*30), ...
        'horizontalalignment', 'center', 'handlevisibility', 'off');

    if i == max(size(th))
        loc = int2str(0);
    else
        loc = int2str(180 + i*30);
    end

    %--- Write text in the opposite half of the plot ----------------------
    text(-rt*snt(i), -rt*cst(i), loc, ...
        'handlevisibility', 'off', 'horizontalalignment', 'center');
end

%% Plot elevation grid ====================================================

%--- Define a "unit" radius circle ----------------------------------------
th = 0 : pi/50 : 2*pi;
xunit = cos(th);
yunit = sin(th);

%--- Plot elevation grid lines and tick text ------------------------------
for elevation = 0 : 15 : 90
    elevationSpherical = 90*cos((pi/180) * elevation);

    line(yunit * elevationSpherical, xunit * elevationSpherical, ...
        'lineStyle', ':', 'color', tc, 'linewidth', 0.5, ...
        'handlevisibility', 'off');

    text(0, elevationSpherical, num2str(elevation), ...
        'BackgroundColor', 'white', 'horizontalalignment','center', ...
        'handlevisibility', 'off');
end

%--- Set view to 2-D ------------------------------------------------------
view(0, 90);

%--- Set axis limits ------------------------------------------------------
%save some space for the title
axis([-95 95 -90 101]) ;

for s=1:32
    % for each sattelite
    satAz = Az(s,sMap(s,:)~=0)/180*pi ;
    satEl = El(s,sMap(s,:)~=0)/180*pi ;
    if ~isempty(satAz)

        % Transform elevation angle to a distance to the center of the plot ------
        elSpherical = 90*cos(satEl) ;
        
        %--- Transform data to Cartesian coordinates ------------------------------
        yy = elSpherical .* cos( satAz ) ;
        xx = elSpherical .* sin( satAz ) ;
        
        plot(hAxis, xx', yy','LineWidth',2,'Color',[0.2 0.4 0.9]) ;
        
        plot(hAxis, xx(end), yy(end), 'rs', 'MarkerSize', 15) ;

        text(xx(end), yy(end), sprintf('%1d',s), 'HorizontalAlignment','Center','FontSize',12) ;
        
    end
end

%--- Make sure both axis have the same data aspect ratio ------------------
axis(hAxis, 'equal') ;

%--- Switch off the standard Cartesian axis -------------------------------
axis(hAxis, 'off');