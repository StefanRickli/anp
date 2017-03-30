function [] = d_shape_06a_1st_interval()
% do the first element
    % can be
    % either: a negative imaginary pz with overlap up to a positive one that
    %         just touches the origin
    % or: the first straight interval
    
    interval_ii = 1;
    idx_current_pz = NaN;
    pz = [im_pz_sorted.value];
    idx_first_positive = find(pz >= 0,1,'first');
    idx_last_negative = find(pz < 0,1,'last');
    positive_pz_remain = length(pz)-idx_first_positive+1;
    prev_upper_bound = 0;
    
    if any([im_pz_sorted.neg_overlapping])
            idx_current_pz = idx_last_negative;
            
            interval_list(1).type = get_detour_type(im_pz_sorted(idx_current_pz).type,'part');
            interval_list(1).q(1) = 0;
            
            arc_length_overlapping_pole = radii.detour_pole * (asin(-im_pz_sorted(idx_current_pz).value/radii.detour_pole) + angles.detour_pole_phi0);
            arc_length_overlapping_zero = radii.detour_zero * (asin(-im_pz_sorted(idx_current_pz).value/radii.detour_zero) + angles.detour_zero_phi0);
            interval_length = im_pz_sorted(idx_current_pz).pole*(arc_lengths.detour_pole - arc_length_overlapping_pole) + im_pz_sorted(idx_current_pz).zero*(arc_lengths.detour_zero - arc_length_overlapping_zero);
            interval_list(1).q_len = interval_length;
            
            interval_list(1).q(2) = interval_length;
            prev_upper_bound = interval_list(1).q(2);
            
            switch im_pz_sorted(end).type
                case 'p'
                    interval_list(1).input_fct_handle = @(q) circ_detour(map(q,0,interval_list(1).q(2),arc_lengths.detour_pole-interval_list(1).q_len,arc_lengths.detour_pole),radii.detour_pole,secant_pole,im_pz_sorted(idx_current_pz).value);
                case 'z'
                    interval_list(1).input_fct_handle = @(q) circ_detour(map(q,0,interval_list(1).q(2),arc_lengths.detour_zero-interval_list(1).q_len,arc_lengths.detour_zero),radii.detour_zero,secant_zero,im_pz_sorted(idx_current_pz).value);
            end
            
            
            interval_ii = interval_ii + 1;
            
            
        tools.dbg('negative and overlapping p/z\n');
        tools.dbg('interval\t[%.3f\t%.3f],\tlength = %.3f,\tdetour\n',interval_list(1).q(1),interval_list(1).q(2),interval_list(1).q(2)-interval_list(1).q(1));
        
        
    elseif any([im_pz_sorted.pos_overlapping])
            idx_current_pz = idx_first_positive;
            
            interval_list(1).type = get_detour_type(im_pz_sorted(idx_current_pz).type,'part');
            interval_list(1).q(1) = 0;
            
            arc_length_overlapping_pole = radii.detour_pole * (asin(-im_pz_sorted(idx_current_pz).value/radii.detour_pole) + angles.detour_pole_phi0);
            arc_length_overlapping_zero = radii.detour_zero * (asin(-im_pz_sorted(idx_current_pz).value/radii.detour_zero) + angles.detour_zero_phi0);
            interval_length = im_pz_sorted(idx_current_pz).pole*(arc_lengths.detour_pole - arc_length_overlapping_pole) + im_pz_sorted(idx_current_pz).zero*(arc_lengths.detour_zero - arc_length_overlapping_zero);
            interval_list(1).q_len = interval_length;
            
            interval_list(1).q(2) = interval_length;
            prev_upper_bound = interval_list(1).q(2);
            
            switch im_pz_sorted(idx_current_pz).type
                case 'p'
                    interval_list(1).input_fct_handle = @(q) circ_detour(map(q,0,interval_list(1).q(2),arc_lengths.detour_pole-interval_list(1).q_len,arc_lengths.detour_pole),radii.detour_pole,secant_pole,im_pz_sorted(idx_current_pz).value);
                case 'z'
                    interval_list(1).input_fct_handle = @(q) circ_detour(map(q,0,interval_list(1).q(2),arc_lengths.detour_zero-interval_list(1).q_len,arc_lengths.detour_zero),radii.detour_zero,secant_zero,im_pz_sorted(idx_current_pz).value);
            end

            positive_pz_remain = positive_pz_remain - 1;
            
            interval_ii = interval_ii + 1;
            
            
        tools.dbg('positive and overlapping p/z\n');
        tools.dbg('interval\t[%.3f\t%.3f],\tlength = %.3f,\tdetour\n',interval_list(1).q(1),interval_list(1).q(2),interval_length);
        
        
    elseif any([im_pz_sorted.pos_on_origin])
            idx_current_pz = idx_first_positive;
            
            interval_list(1).type = get_detour_type(im_pz_sorted(idx_current_pz).type,[]);
            interval_list(1).q(1) = 0;
            
            interval_length = im_pz_sorted(idx_current_pz).pole*arc_lengths.detour_pole + im_pz_sorted(idx_current_pz).zero*arc_lengths.detour_zero;
            interval_list(1).q_len = interval_length;
            
            interval_list(1).q(2) = interval_length;
            prev_upper_bound = interval_list(1).q(2);
            
            switch im_pz_sorted(idx_current_pz).type
                case 'p'
                    interval_list(interval_ii).input_fct_handle = @(q) circ_detour(map(q,interval_list(interval_ii).q(1),interval_list(interval_ii).q(2),0,arc_lengths.detour_pole),radii.detour_pole,secant_pole,im_pz_sorted(idx_current_pz).value);
                case 'z'
                    interval_list(interval_ii).input_fct_handle = @(q) circ_detour(map(q,interval_list(interval_ii).q(1),interval_list(interval_ii).q(2),0,arc_lengths.detour_zero),radii.detour_zero,secant_zero,im_pz_sorted(idx_current_pz).value);
            end
            
            positive_pz_remain = positive_pz_remain - 1;
            
            interval_ii = interval_ii + 1;
            
            
        tools.dbg('positive p/z that ends right on origin\n');
        tools.dbg('interval\t[%.3f\t%.3f],\tlength = %.3f,\tlinear\n',interval_list(1).q(1),interval_list(1).q(2),interval_length);
        
        
    elseif  ~isempty(idx_first_positive)
            prev_upper_bound = 0;
%             interval_list(1).q(1) = 0;
%             
%             interval_list(1).q_len = im_pz_sorted(positive_idx).value - halfsecant_pole;
%             
%             interval_list(1).q(2) = im_pz_sorted(positive_idx).value - halfsecant_pole;
%             prev_upper_bound = interval_list(1).q(2);
        tools.dbg('no p/z in range of origin\n');
        
        
    else
            prev_upper_bound = 0;
        
        tools.dbg('no imag positive p/z\n');
        
        
    end
    
    tools.dbg('------------------------------------------------------------------\n');
end