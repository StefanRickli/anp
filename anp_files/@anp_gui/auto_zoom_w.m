% The purpose of this function is to provide a means of detecting clumping
% of zeros and poles. If we have a bunch of zeros and poles that are fairly
% close to each other with one outlier, forget about the outlier and zoom
% in on the group.
function [w_xlim,w_ylim] = anp_plot_auto_zoom_w(w_values)
    % for the right subplot it's fairly easy. Have the Plot such that the
    % whole nyquist curve fits into it.
    w_xlim = anp_stretch_centered([min([real(w_values),0]),max([real(w_values),0])],1.15);
    w_ylim = anp_stretch_centered([min([imag(w_values),0]),max([imag(w_values),0])],1.15);
    
    width = w_xlim(2) - w_xlim(1);
    height = w_ylim(2) - w_ylim(1);

    % make this plot also rectangular
    if width >= height
        y0 = mean(w_ylim);
        w_ylim = [(y0-width/2),(y0+width/2)];
    else
        x0 = mean(w_xlim);
        w_xlim = [(x0-height/2),(x0+height/2)];
    end
end
    