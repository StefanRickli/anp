function [interval_ii,idx_current_pz,prev_upper_bound] = ds07g_first_after_crop2(this,interval_ii)
    im_pz_sorted =          this.im_pz_sorted;
    interval_list =         this.interval_list;
    positions =             this.positions;
    halfsecant_pole =       this.halfsecant_pole;
    halfsecant_zero =       this.halfsecant_zero;
    
    % do the first linear interval after crop2
    interval_list(interval_ii).type = 'axis';
    interval_list(interval_ii).q(1) = interval_list(interval_ii-1).q(2);
    
    if any([im_pz_sorted.value] < 0)
        idx_current_pz = 1;
        interval_length = positions.crop_y0 - (-im_pz_sorted(idx_current_pz).value + im_pz_sorted(idx_current_pz).pole*halfsecant_pole + im_pz_sorted(idx_current_pz).zero*halfsecant_zero);
        zb = im_pz_sorted(idx_current_pz).value - im_pz_sorted(idx_current_pz).pole*halfsecant_pole - im_pz_sorted(idx_current_pz).zero*halfsecant_zero;
    elseif any([im_pz_sorted.pos_overlapping])
        idx_current_pz = NaN;
        interval_length = positions.crop_y0 - (-im_pz_sorted(1).value + im_pz_sorted(1).pole*halfsecant_pole + im_pz_sorted(1).zero*halfsecant_zero);
        zb = im_pz_sorted(1).value - im_pz_sorted(1).pole*halfsecant_pole - im_pz_sorted(1).zero*halfsecant_zero;
    else
        idx_current_pz = NaN;
        interval_length = positions.crop_y0;
        zb = 0;
    end
    interval_list(interval_ii).q_len = interval_length;
    
    interval_list(interval_ii).q(2) = interval_list(interval_ii).q(1) + interval_length;
    prev_upper_bound = interval_list(interval_ii).q(2);
    
    za = -positions.crop_y0;
    interval_list(interval_ii).input_fct_handle = @(q) im_axis_line(q,interval_list(interval_ii).q(1),interval_list(interval_ii).q(2),za,zb);
    
    tools.dbg('interval\t[%.3f\t%.3f],\tlength = %.3f,\tlinear_first_pos_after_crop2\n',interval_list(interval_ii).q(1),interval_list(interval_ii).q(2),interval_length);
    interval_ii = interval_ii + 1;
    
    this.interval_list = interval_list;
end