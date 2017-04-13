% The purpose of this function is to provide a means of detecting clumping
% of zeros and poles. If we have a bunch of zeros and poles that are fairly
% close to each other with one outlier, forget about the outlier and zoom
% in on the group.
function [z_xlim,z_ylim] = anp_plot_auto_zoom_z(zp, R)
    % find the 'center of gravity' of the zeros and poles
    locus_m = mean(zp);
    
    % where are 66% of the points?
    locus_std_real = std(real(zp));
    locus_std_imag = std(imag(zp));
    
    % get indexes of those points that lie within the ellipse given by the stddev in x- and y-direction
    locus_inner =  find(abs(zp-locus_m) <= sqrt(locus_std_real^2+locus_std_imag^2));
    
    % how much area do these points cover in the zoomed-out plot (when the D-curve fits into it)?
    locus_inner_dims = [max(real(zp(locus_inner)))-min(real(zp(locus_inner)));max(imag(zp(locus_inner)))-min(imag(zp(locus_inner)))];
    if min(locus_inner_dims)/max(locus_inner_dims) <= 1/50
        % almost one-dimensional
        locus_inner_dim = max(locus_inner_dims);
        locus_fract = locus_inner_dim/(2.01*R);
    else
        locus_inner_area = locus_inner_dims(1)*locus_inner_dims(2);
        locus_fract = locus_inner_area/(2.01*R)^2;
    end
    
    % if that area is below 0.2% then we zoom into the bulk of the points,
    % leaving out outliers
    if locus_fract > 0.015
        % all good, fit the whole D-curve and the zp-points into the plot
        z_xlim = 1.05*R*[-1,1];
        z_ylim = 1.05*R*[-1,1];
    else
        % include the origin in the box around the clumped points and then
        % make the plot size 8 times (rule of thumb) larger
        z_xlim = anp_stretch_centered([min([real(zp(locus_inner)),0]),max([real(zp(locus_inner)),0])],8);
        z_ylim = anp_stretch_centered([min([imag(zp(locus_inner)),0]),max([imag(zp(locus_inner)),0])],8);
        
        width = z_xlim(2) - z_xlim(1);
        height = z_ylim(2) - z_ylim(1);
        
        % make the plot rectangular (make the smaller axis as large as the
        % bigger one)
        if width >= height
            y0 = mean(z_ylim);
            z_ylim = [(y0-width/2),(y0+width/2)];
        else
            x0 = mean(z_xlim);
            z_xlim = [(x0-height/2),(x0+height/2)];
        end
    end
end
    