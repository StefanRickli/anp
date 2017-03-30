function [] = d_shape_06h_negative_pz()
    
        % treat negative p/z
    negative_pz_remain = sum(pz<0) - any([[im_pz_sorted.neg_on_origin],[im_pz_sorted.neg_overlapping]]);
    
    while(negative_pz_remain)
        % there are negative poles/zeros left to treat
        
        interval_list(interval_ii).type = get_detour_type(im_pz_sorted(idx_current_pz).type,[]);
        interval_list(interval_ii).q(1) = prev_upper_bound;
        
        interval_length = im_pz_sorted(idx_current_pz).pole*arc_lengths.detour_pole + im_pz_sorted(idx_current_pz).zero*arc_lengths.detour_zero;
        interval_list(interval_ii).q_len = interval_length;
        
        interval_list(interval_ii).q(2) = interval_list(interval_ii).q(1) + interval_length;
        prev_upper_bound = interval_list(interval_ii).q(2);
        
        switch im_pz_sorted(idx_current_pz).type
            case 'p'
                interval_list(interval_ii).input_fct_handle = @(q) circ_detour(map(q,interval_list(interval_ii).q(1),interval_list(interval_ii).q(2),0,arc_lengths.detour_pole),radii.detour_pole,secant_pole,im_pz_sorted(idx_current_pz).value);
            case 'z'
                interval_list(interval_ii).input_fct_handle = @(q) circ_detour(map(q,interval_list(interval_ii).q(1),interval_list(interval_ii).q(2),0,arc_lengths.detour_zero),radii.detour_zero,secant_zero,im_pz_sorted(idx_current_pz).value);
        end
        
        tools.dbg('interval\t[%.3f\t%.3f],\tlength = %.3f,\tdetour_neg\n',interval_list(interval_ii).q(1),interval_list(interval_ii).q(2),interval_length);
        interval_ii = interval_ii + 1;
        
        if negative_pz_remain > 1
            interval_list(interval_ii).type = 'axis';
            interval_list(interval_ii).q(1) = prev_upper_bound;
            
            interval_length = pole_zero_combinations(idx_current_pz).distance - halfsecant_pole*(im_pz_sorted(idx_current_pz).pole + im_pz_sorted(idx_current_pz+1).pole) - halfsecant_zero*(im_pz_sorted(idx_current_pz).zero + im_pz_sorted(idx_current_pz+1).zero);            
            interval_list(interval_ii).q_len = interval_length;
            
            interval_list(interval_ii).q(2) = interval_list(interval_ii).q(1) + interval_length;
            prev_upper_bound = interval_list(interval_ii).q(2);
            
            za = im_pz_sorted(idx_current_pz).value + im_pz_sorted(idx_current_pz).pole*halfsecant_pole + im_pz_sorted(idx_current_pz).zero*halfsecant_zero;
            zb = im_pz_sorted(idx_current_pz+1).value - im_pz_sorted(idx_current_pz+1).pole*halfsecant_pole - im_pz_sorted(idx_current_pz+1).zero*halfsecant_zero;
            interval_list(interval_ii).input_fct_handle = @(q) im_axis_line(q,interval_list(interval_ii).q(1),interval_list(interval_ii).q(2),za,zb);

            tools.dbg('interval\t[%.3f\t%.3f],\tlength = %.3f,\tlinear_neg\n',interval_list(interval_ii).q(1),interval_list(interval_ii).q(2),interval_length);
            idx_current_pz = idx_current_pz + 1;
            interval_ii = interval_ii + 1;
            negative_pz_remain = negative_pz_remain - 1;
        else
            break;
        end
    end
end