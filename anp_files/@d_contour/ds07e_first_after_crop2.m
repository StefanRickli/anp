function [interval_ii,idx_current_pz,prev_upper_bound] = ds07e_first_after_crop2(this,interval_ii)
    % Handles the first interval after the second crop circle.
    % This is certainly a piece of straight Im-axis.
    
    % Make the used variables local to get rid of 'this.'
    im_pz_sorted =          this.im_pz_sorted;
    interval_list =         this.interval_list;
    positions =             this.positions;
    halfsecant_pole =       this.halfsecant_pole;
    halfsecant_zero =       this.halfsecant_zero;
    
    % Note down the type of this interval.
    interval_list(interval_ii).type = 'axis';
    
    % Note down where this interval starts on the whole length of
    % the D-contour.
    interval_list(interval_ii).q(1) = interval_list(interval_ii-1).q(2);
    
    % At this point there are three possibilities up to where the current
    % interval goes. They will be encountered in the order of the following
    % IF-statements:
    if any([im_pz_sorted.value] < 0)
        % Hit the detour of a p/z that lies on the lower HP of the Im-axis.
        
        % As the list of p/z is sorted by increasing imaginary part, we
        % take the first to appear on the list.
        idx_current_pz = 1;
        
        % The length of this interval is the distance from the end of
        % Crop 2 up to the beginning of this first detour.
        interval_length = positions.crop_y0 + ...
                          -( -im_pz_sorted(idx_current_pz).value + ...
                             halfsecant_pole * im_pz_sorted(idx_current_pz).pole + ...
                             halfsecant_zero * im_pz_sorted(idx_current_pz).zero ...
                            );
        
        % 'zb' marks the end point where the straight axis part stops.
        % Here it's where the first detour starts.
        zb = im_pz_sorted(idx_current_pz).value + ...
             -halfsecant_pole * im_pz_sorted(idx_current_pz).pole + ...
             -halfsecant_zero * im_pz_sorted(idx_current_pz).zero;
         
    elseif any([im_pz_sorted.pos_overlapping])
        % An imaginary p/z of the upper HP that lies so close to the origin
        % that its detour extends beyond the real axis into the lower HP.
        
        % Signal to the next function that no negative p/z 
        idx_current_pz = NaN;
        interval_length = positions.crop_y0 + ...
                          -( -im_pz_sorted(1).value + ...
                             halfsecant_pole * im_pz_sorted(1).pole + ...
                             halfsecant_zero * im_pz_sorted(1).zero ...
                            );
        zb = im_pz_sorted(1).value + ...
             -halfsecant_pole * im_pz_sorted(1).pole + ...
             -halfsecant_zero * im_pz_sorted(1).zero;
    else
        % None of the above cases: 
        
        idx_current_pz = NaN;
        interval_length = positions.crop_y0;
        zb = 0;
    end
    interval_list(interval_ii).q_len = interval_length;
    
    interval_list(interval_ii).q(2) = interval_list(interval_ii).q(1) + interval_length;
    prev_upper_bound = interval_list(interval_ii).q(2);
    
    za = -positions.crop_y0;
    interval_list(interval_ii).input_fct_handle = @(q) im_axis_line(q,interval_list(interval_ii).q(1),interval_list(interval_ii).q(2),za,zb);
    
    dbg_out('interval\t[%.3f\t%.3f],\tlength = %.3f,\tlinear_first_pos_after_crop2\n',interval_list(interval_ii).q(1),interval_list(interval_ii).q(2),interval_length);
    interval_ii = interval_ii + 1;
    
    this.interval_list = interval_list;
end