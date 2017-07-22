% The purpose of this function is to provide a means of detecting clumping
% of zeros and poles. If we have a bunch of zeros and poles that are fairly
% close to each other with one outlier, forget about the outlier and zoom
% in on the group.
function [w_xlim,w_ylim] = auto_zoom_w(this, w_values)
    % for the right subplot it's fairly easy. Have the Plot such that the
    % whole nyquist curve fits into it.
    w_xlim = this.stretch_centered([min([real(w_values),0]),max([real(w_values),0])],1.15);
    w_ylim = this.stretch_centered([min([imag(w_values),0]),max([imag(w_values),0])],1.15);
end
    