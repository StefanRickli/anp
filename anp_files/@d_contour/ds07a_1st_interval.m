function [interval_ii,idx_current_pz,prev_upper_bound,positive_pz_remain,idx_first_positive,idx_last_negative] = ...
        ds07a_1st_interval(this)
    % Handles the first interval in the list.
    % This could be: - a regular case: a straight line part
    %                                  --> Let '..._positive_pz.m' handle
    %                                      this.
    %                - the upper part of a detour that crosses the real
    %                  axis.
    %                - the detour of a p/z that exactly ends on the origin
    
    % Make the used variables local to get rid of 'this.'
    im_pz_sorted =          this.im_pz_sorted;
    interval_list =         this.interval_list;
    radii =                 this.radii;
    arc_lengths =           this.arc_lengths;
    secant_pole =           this.secant_pole;
    secant_zero =           this.secant_zero;
    angles =                this.angles;
    
    % Init iterables
    interval_ii =           1;
    idx_current_pz =        NaN;
    prev_upper_bound =      0;
    
    % Get indexes of the first positive and last negative purely imaginary
    % pole / zero (if any exist).
    pz = [im_pz_sorted.value];
    idx_first_positive =    find(pz >= 0,1,'first');
    idx_last_negative =     find(pz < 0,1,'last');
    
    % Used for loop condition later in the code.
    positive_pz_remain =    length(pz)-idx_first_positive+1;
    
    if any([im_pz_sorted.neg_overlapping])
        % Case 1/5: There is a negative p/z whose detour crosses the real
        %           axis towards the upper HP.
        
            idx_current_pz = idx_last_negative;
            
            % Note down the type of this interval.
            interval_list(1).type = [repmat('detour_pole_part', im_pz_sorted(idx_current_pz).pole), ...
                                     repmat('detour_zero_part', im_pz_sorted(idx_current_pz).zero)];
            
            % Note down where this interval starts on the whole length of
            % the D-contour.
            interval_list(1).q(1) = 0;
            
            % Note down the effective length this interval takes.
            % If it's a straight line, then q_len is just |begin - end|,
            % if it's part of an arc, then q_len is the length of that arc.
            arc_length_overlapping_pole = radii.detour_pole * ( asin(-im_pz_sorted(idx_current_pz).value/radii.detour_pole) + ...
                                                                angles.detour_pole_phi0 ...
                                                               );
            arc_length_overlapping_zero = radii.detour_zero * ( asin(-im_pz_sorted(idx_current_pz).value/radii.detour_zero) + ...
                                                                angles.detour_zero_phi0 ...
                                                               );
            % (Remember that im_pz_sorted(ii).pole and .zero are booleans)
            interval_length = (arc_lengths.detour_pole - arc_length_overlapping_pole) * im_pz_sorted(idx_current_pz).pole + ...
                              (arc_lengths.detour_zero - arc_length_overlapping_zero) * im_pz_sorted(idx_current_pz).zero;
            interval_list(1).q_len = interval_length;
            
            % Note down where this interval ends on the whole length of
            % the D-contour.
            interval_list(1).q(2) = interval_length;
            
            % Remember the end of this interval for ease of use later.
            prev_upper_bound = interval_list(1).q(2);
            
            % Construct a function handle which will be fed with q = 0 to
            % q = 'length of this interval'.
            switch im_pz_sorted(idx_current_pz).type
                case 'p'
                    interval_list(1).input_fct_handle = @(q) circ_detour(map(q,0,interval_list(1).q(2),arc_lengths.detour_pole-interval_list(1).q_len,arc_lengths.detour_pole), ...
                                                                         radii.detour_pole, ...
                                                                         secant_pole, ...
                                                                         im_pz_sorted(idx_current_pz).value);
                case 'z'
                    interval_list(1).input_fct_handle = @(q) circ_detour(map(q,0,interval_list(1).q(2),arc_lengths.detour_zero-interval_list(1).q_len,arc_lengths.detour_zero), ...
                                                                         radii.detour_zero, ...
                                                                         secant_zero, ...
                                                                         im_pz_sorted(idx_current_pz).value);
            end
            
            
            interval_ii = interval_ii + 1;
            
            
        dbg_out('negative and overlapping p/z\n');
        dbg_out('interval\t[%.3f\t%.3f],\tlength = %.3f,\tdetour\n',interval_list(1).q(1),interval_list(1).q(2),interval_list(1).q(2)-interval_list(1).q(1));
        
        
    elseif any([im_pz_sorted.pos_overlapping])
        % Case 2/5: There is a positive p/z whose detour crosses the real
        %           axis towards the lower HP.
        
            idx_current_pz = idx_first_positive;
            
            % Note down the type of this interval.
            interval_list(1).type = [repmat('detour_pole_part', im_pz_sorted(idx_current_pz).pole), ...
                                     repmat('detour_zero_part', im_pz_sorted(idx_current_pz).zero)];
            
            % Note down where this interval starts on the whole length of
            % the D-contour.
            interval_list(1).q(1) = 0;
            
            % Note down the effective length this interval takes.
            % If it's a straight line, then q_len is just |begin - end|,
            % if it's part of an arc, then q_len is the length of that arc.
            arc_length_overlapping_pole = radii.detour_pole * ( asin(-im_pz_sorted(idx_current_pz).value/radii.detour_pole) + ...
                                                                angles.detour_pole_phi0 ...
                                                               );
            arc_length_overlapping_zero = radii.detour_zero * ( asin(-im_pz_sorted(idx_current_pz).value/radii.detour_zero) + ...
                                                                angles.detour_zero_phi0 ...
                                                               );
            % (Remember that im_pz_sorted(ii).pole and .zero are booleans)
            interval_length = (arc_lengths.detour_pole - arc_length_overlapping_pole) * im_pz_sorted(idx_current_pz).pole + ...
                              (arc_lengths.detour_zero - arc_length_overlapping_zero) * im_pz_sorted(idx_current_pz).zero;
            interval_list(1).q_len = interval_length;
            
            % Note down where this interval ends on the whole length of
            % the D-contour.
            interval_list(1).q(2) = interval_length;
            
            % Remember the end of this interval for ease of use later.
            prev_upper_bound = interval_list(1).q(2);
            
            % Construct a function handle which will be fed with 0 to
            % 'length of this interval'.
            switch im_pz_sorted(idx_current_pz).type
                case 'p'
                    interval_list(1).input_fct_handle = @(q) circ_detour(map(q,0,interval_list(1).q(2),arc_lengths.detour_pole-interval_list(1).q_len,arc_lengths.detour_pole), ...
                                                                         radii.detour_pole, ...
                                                                         secant_pole, ...
                                                                         im_pz_sorted(idx_current_pz).value);
                case 'z'
                    interval_list(1).input_fct_handle = @(q) circ_detour(map(q,0,interval_list(1).q(2),arc_lengths.detour_zero-interval_list(1).q_len,arc_lengths.detour_zero), ...
                                                                         radii.detour_zero, ...
                                                                         secant_zero, ...
                                                                         im_pz_sorted(idx_current_pz).value);
            end

            positive_pz_remain = positive_pz_remain - 1;
            
            interval_ii = interval_ii + 1;
            
            
        dbg_out('positive and overlapping p/z\n');
        dbg_out('interval\t[%.3f\t%.3f],\tlength = %.3f,\tdetour\n',interval_list(1).q(1),interval_list(1).q(2),interval_length);
        
        
    elseif any([im_pz_sorted.pos_on_origin])
        % Case 3/5: There is a positive p/z whose detour ends exactly on
        %           the origin.
        
            idx_current_pz = idx_first_positive;
            
            % Note down the type of this interval.
            interval_list(1).type = [repmat('detour_pole', im_pz_sorted(idx_current_pz).pole), ...
                                     repmat('detour_zero', im_pz_sorted(idx_current_pz).zero)];
            
            % Note down where this interval starts on the whole length of
            % the D-contour.
            interval_list(1).q(1) = 0;
            
            % Note down the effective length this interval takes.
            % If it's a straight line, then q_len is just |begin - end|,
            % if it's part of an arc, then q_len is the length of that arc.
            interval_length = arc_lengths.detour_pole * im_pz_sorted(idx_current_pz).pole + ...
                              arc_lengths.detour_zero * im_pz_sorted(idx_current_pz).zero;
            interval_list(1).q_len = interval_length;
            
            % Note down where this interval ends on the whole length of
            % the D-contour.
            interval_list(1).q(2) = interval_length;
            
            % Remember the end of this interval for ease of use later.
            prev_upper_bound = interval_list(1).q(2);
            
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
            
            positive_pz_remain = positive_pz_remain - 1;
            
            interval_ii = interval_ii + 1;
            
            
        dbg_out('positive p/z that ends right on origin\n');
        dbg_out('interval\t[%.3f\t%.3f],\tlength = %.3f,\tlinear\n',interval_list(1).q(1),interval_list(1).q(2),interval_length);
        
        
    elseif  ~isempty(idx_first_positive)
        % Case 4/5: Although there's at least one purely imaginary p/z on
        %           the positive HP, its detour doesn't get near the real
        %           axis.
        
            prev_upper_bound = 0;
            
        dbg_out('no p/z in range of origin\n');
        
        
    else
        % Case 5/5: There are no purely imaginary p/z.
        
            prev_upper_bound = 0;
        
        dbg_out('no imag positive p/z\n');
        
        
    end
    
    dbg_out('------------------------------------------------------------------\n');
    
    % Write back changes to the interval list.
    this.interval_list = interval_list;
end