function [interval_ii] = ds07c_last_straight_before_crop1(this,interval_ii,idx_current_pz,prev_upper_bound)
    im_pz_sorted =          this.im_pz_sorted;
    interval_list =         this.interval_list;
    halfsecant_pole =       this.halfsecant_pole;
    halfsecant_zero =       this.halfsecant_zero;
    positions =             this.positions;
    
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
    
    dbg_out('interval\t[%.3f\t%.3f],\tlength = %.3f,\tlinear_last_pos_before_crop1\n',interval_list(interval_ii).q(1),interval_list(interval_ii).q(2),interval_length);
    interval_ii = interval_ii + 1;
    
    this.interval_list = interval_list;
end