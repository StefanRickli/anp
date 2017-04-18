function [interval_ii,idx_current_pz,prev_upper_bound] = ds07b_positive_pz(this,interval_ii,idx_current_pz,prev_upper_bound,positive_pz_remain,idx_first_positive)
    % Handles the intervals that deal with the positive imaginary axis up until the last straight part.
    % The shape functions can be: - straight line parts
    %                             - detour parts
    
    % Make the used variables local to get rid of 'this.'
    im_pz_sorted =          this.im_pz_sorted;
    pole_zero_combinations =this.im_pz_combinations;
    interval_list =         this.interval_list;
    secant_pole =           this.secant_pole;
    secant_zero =           this.secant_zero;
    halfsecant_pole =       this.halfsecant_pole;
    halfsecant_zero =       this.halfsecant_zero;
    arc_lengths =           this.arc_lengths;
    radii =                 this.radii;
    
    % Every imaginary p/z on the axis at this point is preceded by a piece
    % of straight line on the Im-axis.
    while(positive_pz_remain)
        % There are pure imaginary p/z on the positive HP left to treat
        
        idx_current_pz = idx_current_pz + 1;
        
        % Step 1: Define a piece of straight Im-axis part up to the detour
        %         arc of the next p/z
        
        
        % Note down the type of this interval.
        interval_list(interval_ii).type = 'axis';
        
        % Note down where this interval starts on the whole length of
        % the D-contour.
        interval_list(interval_ii).q(1) = prev_upper_bound;
        
        if isnan(idx_current_pz)
            % 'idx_current_pz' is NaN when '..._1st_interval.m' didn't do
            % anything. That is the case when no detour of the positive HP
            % lands on the origin or no detour crosses the real axis.
            
            % Init the iterator 'idx_current_pz'.
            idx_current_pz = idx_first_positive;
            
            % The straight axis part is simply the distance from the origin
            % to the first detour arc beginning.
            interval_length = im_pz_sorted(idx_current_pz).value - ...
                              halfsecant_pole * im_pz_sorted(idx_current_pz).pole - ...
                              halfsecant_zero * im_pz_sorted(idx_current_pz).zero;
            
            % 'za' marks the start point of the straight axis part.
            % Here it's at the origin.
            za = 0;
        else
            % 1st entry of the interval list was a special case, i.e. we
            % began with (part of) a detour. So the length of the first
            % straight axis part is not from origin but the end of the
            % detour.
            
            interval_length = pole_zero_combinations(idx_current_pz-1).distance - ...
                              halfsecant_pole * ( im_pz_sorted(idx_current_pz-1).pole + im_pz_sorted(idx_current_pz).pole ) - ...
                              halfsecant_zero * ( im_pz_sorted(idx_current_pz-1).zero + im_pz_sorted(idx_current_pz).zero );            
            
            % 'za' marks the start point of the straight axis part.
            % Here it's at the end of the preceeding detour.
            za = im_pz_sorted(idx_current_pz-1).value + ...
                 halfsecant_pole * im_pz_sorted(idx_current_pz-1).pole + ...
                 halfsecant_zero * im_pz_sorted(idx_current_pz-1).zero;
        end
        
        % Note down the length this interval takes.
        interval_list(interval_ii).q_len = interval_length;
        
        % Note down where this interval ends on the whole length of
        % the D-contour.
        interval_list(interval_ii).q(2) = interval_list(interval_ii).q(1) + interval_length;
        
        % Remember the end of this interval for ease of use later.
        prev_upper_bound = interval_list(interval_ii).q(2);
        
        % 'zb' marks the end point where the straight axis part stops.
        % Here it's where the next detour starts.
        zb = im_pz_sorted(idx_current_pz).value - ...
             halfsecant_pole * im_pz_sorted(idx_current_pz).pole - ...
             halfsecant_zero * im_pz_sorted(idx_current_pz).zero;
        
        % Construct a function handle which will be fed with q = 0 to
        % q = 'length of this interval'.
        interval_list(interval_ii).input_fct_handle = @(q) im_axis_line(q,interval_list(interval_ii).q(1),interval_list(interval_ii).q(2),za,zb);
        
        dbg_out('interval\t[%.3f\t%.3f],\tlength = %.3f,\tlinear_pos\n',interval_list(interval_ii).q(1),interval_list(interval_ii).q(2),interval_length);
        interval_ii = interval_ii + 1;
        
        % -----------------------------------------------------------------
        
        % Step 2: Define the detour after the piece of straight
        %         Im-axis part.
        
        % Note down the type of this interval.
        interval_list(interval_ii).type = [repmat('detour_pole', im_pz_sorted(idx_current_pz).pole), ...
                                           repmat('detour_zero', im_pz_sorted(idx_current_pz).zero)];
        
        % Note down where this interval starts on the whole length of
        % the D-contour.
        interval_list(interval_ii).q(1) = prev_upper_bound;
        
        % Note down the length this interval takes. In this case it's part
        % of a detour arc, so q_len is the length of that arc.
        interval_length = arc_lengths.detour_pole * im_pz_sorted(idx_current_pz).pole + ...
                          arc_lengths.detour_zero * im_pz_sorted(idx_current_pz).zero;
        interval_list(interval_ii).q_len = interval_length;
        
        % Note down where this interval ends on the whole length of
        % the D-contour.
        interval_list(interval_ii).q(2) = interval_list(interval_ii).q(1) + interval_length;
        
        % Remember the end of this interval for ease of use later.
        prev_upper_bound = interval_list(interval_ii).q(2);
        
        % Construct a function handle which will be fed with q = 0 to
        % q = 'length of this interval'.
        switch im_pz_sorted(idx_current_pz).type
            case 'p'
                interval_list(interval_ii).input_fct_handle = @(q) circ_detour(map(q,interval_list(interval_ii).q(1),interval_list(interval_ii).q(2),0,arc_lengths.detour_pole), ...
                                                                               radii.detour_pole, ...
                                                                               secant_pole, ...
                                                                               im_pz_sorted(idx_current_pz).value);
            case 'z'
                interval_list(interval_ii).input_fct_handle = @(q) circ_detour(map(q,interval_list(interval_ii).q(1),interval_list(interval_ii).q(2),0,arc_lengths.detour_zero), ...
                                                                               radii.detour_zero, ...
                                                                               secant_zero, ...
                                                                               im_pz_sorted(idx_current_pz).value);
        end
        
        dbg_out('interval\t[%.3f\t%.3f],\tlength = %.3f,\tdetour_pos\n',interval_list(interval_ii).q(1),interval_list(interval_ii).q(2),interval_length);
        
        % Update iterator and counting variable.
        interval_ii = interval_ii + 1;
        positive_pz_remain = positive_pz_remain - 1;
    end
    
    % Write back changes to the interval list.
    this.interval_list = interval_list;
end