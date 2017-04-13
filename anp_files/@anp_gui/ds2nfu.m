function [x_t,y_t] = ds2nfu(hAx,x,y,axlim)
    % DS2NFU  Convert data space units into normalized figure units.
    %
    % Michelle Hirsch
    % mhirsch@mathworks.com
    % Copyright 2006-2014 The MathWorks, Inc
    %
    % Shortened by Stefan Rickli
    
    %% Get limits
    axun = get(hAx,'Units');
    set(hAx,'Units','normalized');
    axpos = get(hAx,'Position');
    axwidth = diff(axlim(1:2));
    axheight = diff(axlim(3:4));
    
    %% Transform data
    x_t = (x-axlim(1))*axpos(3)/axwidth + axpos(1);
    y_t = (y-axlim(3))*axpos(4)/axheight + axpos(2);
    
    %% Restore axes units
    set(hAx,'Units',axun)
end
