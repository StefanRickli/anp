function [] = d_shape_06c_last_straight_before_crop1()
    
    % do the last linear interval before crop1
    interval_list(interval_ii).type = 'axis';
    interval_list(interval_ii).q(1) = prev_upper_bound;
    
    if isnan(idx_current_pz)
        interval_length = positions.crop_y0;
        za = 0;
    else
        interval_length = positions.crop_y0 - (im_pz_sorted(idx_current_pz).value + im_pz_sorted(idx_current_pz).pole*halfsecant_pole + im_pz_sorted(idx_current_pz).zero*halfsecant_zero);
        za = im_pz_sorted(idx_current_pz).value + im_pz_sorted(idx_current_pz).pole*halfsecant_pole + im_pz_sorted(idx_current_pz).zero*halfsecant_zero;
    end
    interval_list(interval_ii).q_len = interval_length;
    
    interval_list(interval_ii).q(2) = interval_list(interval_ii).q(1) + interval_length;
    
    zb = positions.crop_y0;
    interval_list(interval_ii).input_fct_handle = @(q) im_axis_line(q,interval_list(interval_ii).q(1),interval_list(interval_ii).q(2),za,zb);
    
    tools.dbg('interval\t[%.3f\t%.3f],\tlength = %.3f,\tlinear_last_pos_before_crop1\n',interval_list(interval_ii).q(1),interval_list(interval_ii).q(2),interval_length);
    interval_ii = interval_ii + 1;
    
    
end