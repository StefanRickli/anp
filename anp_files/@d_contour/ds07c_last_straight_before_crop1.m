function [interval_ii] = ds07c_last_straight_before_crop1(this,interval_ii,idx_current_pz,prev_upper_bound)
    % Handles the last interval before the first crop circle.
    % This is certainly a piece of straight Im-axis.
    
    % Make the used variables local to get rid of 'this.'
    im_pz_sorted =          this.im_pz_sorted;
    interval_list =         this.interval_list;
    halfsecant_pole =       this.halfsecant_pole;
    halfsecant_zero =       this.halfsecant_zero;
    positions =             this.positions;
    
    % Note down the type of this interval.
    interval_list(interval_ii).type = 'axis';
    
    % Note down where this interval starts on the whole length of
    % the D-contour.
    interval_list(interval_ii).q(1) = prev_upper_bound;
    
    if isnan(idx_current_pz)
        % At this point 'idx_current_pz' could still be NaN if there were
        % no detours to handle on the positive HP part of the Im-axis.
        
        % Note down the effective length this interval takes.
        % In this case it's from the origin up to where the upper crop
        % circle starts.
        interval_length = positions.crop_y0;
        
        % 'za' marks the start point of the straight axis part.
        % Here it's the origin.
        za = 0;
    else
        % There was at least one detour on the positive HP part of the
        % Im-axis.
        
        % Note down the effective length this interval takes.
        % In this case it's from where the last detour ended up to where
        % the upper crop circle starts.
        interval_length = positions.crop_y0 + ...
                          -( im_pz_sorted(idx_current_pz).value + ...
                             halfsecant_pole * im_pz_sorted(idx_current_pz).pole + ...
                             halfsecant_zero * im_pz_sorted(idx_current_pz).zero ...
                            );
        
        % 'za' marks the start point of the straight axis part.
        % Here it's where the last detour ended.
        za = im_pz_sorted(idx_current_pz).value + ...
             halfsecant_pole * im_pz_sorted(idx_current_pz).pole + ...
             halfsecant_zero * im_pz_sorted(idx_current_pz).zero;
    end
    
    % Note down the length this interval takes.
    interval_list(interval_ii).q_len = interval_length;
    
    % Note down where this interval ends on the whole length of
    % the D-contour.
    interval_list(interval_ii).q(2) = interval_list(interval_ii).q(1) + interval_length;
    
    % 'zb' marks the end point where the straight axis part stops.
    % Here it's where the upper crop circle starts.
    zb = positions.crop_y0;
    
    % Construct a function handle which will be fed with q = 0 to
    % q = 'length of this interval'.
    interval_list(interval_ii).input_fct_handle = @(q) im_axis_line(q,interval_list(interval_ii).q(1),interval_list(interval_ii).q(2),za,zb);
    
    dbg_out('interval\t[%.3f\t%.3f],\tlength = %.3f,\tlinear_last_pos_before_crop1\n',interval_list(interval_ii).q(1),interval_list(interval_ii).q(2),interval_length);
    
    interval_ii = interval_ii + 1;
    
    % Write back changes to the interval list.
    this.interval_list = interval_list;
end