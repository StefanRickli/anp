function [] = d_shape_06d_crop1_inf_crop2()
    
    % crop1
    interval_list(interval_ii).type = 'crop';
    interval_list(interval_ii).q(1) = interval_list(interval_ii-1).q(2);
    interval_list(interval_ii).q_len = arc_lengths.crop;
    interval_list(interval_ii).q(2) = interval_list(interval_ii).q(1) + arc_lengths.crop;
    
    interval_list(interval_ii).input_fct_handle = @(q) circ_normal(-map(q,interval_list(interval_ii).q(1),interval_list(interval_ii).q(2),0,arc_lengths.crop),radii.crop,pi,positions.crop_x0,positions.crop_y0);
    
    tools.dbg('interval\t[%.3f\t%.3f],\tlength = %.3f,\tcrop1\n',interval_list(interval_ii).q(1),interval_list(interval_ii).q(2),arc_lengths.crop);
    interval_ii = interval_ii + 1;
    
    % inf
    interval_list(interval_ii).type = 'inf';
    interval_list(interval_ii).q(1) = interval_list(interval_ii-1).q(2);
    interval_list(interval_ii).q_len = arc_lengths.inf;
    interval_list(interval_ii).q(2) = interval_list(interval_ii).q(1) + arc_lengths.inf;
    
    interval_list(interval_ii).input_fct_handle = @(q) circ_normal(-map(q,interval_list(interval_ii).q(1),interval_list(interval_ii).q(2),0,arc_lengths.inf),radii.inf,pi/2-angles.crop,0,0);    
    
    tools.dbg('interval\t[%.3f\t%.3f],\tlength = %.3f,\tinf\n',interval_list(interval_ii).q(1),interval_list(interval_ii).q(2),arc_lengths.inf);
    interval_ii = interval_ii + 1;
    
    % crop2
    interval_list(interval_ii).type = 'crop';
    interval_list(interval_ii).q(1) = interval_list(interval_ii-1).q(2);
    interval_list(interval_ii).q_len = arc_lengths.crop;
    interval_list(interval_ii).q(2) = interval_list(interval_ii).q(1) + arc_lengths.crop;
    
    interval_list(interval_ii).input_fct_handle = @(q) circ_normal(-map(q,interval_list(interval_ii).q(1),interval_list(interval_ii).q(2),0,arc_lengths.crop),radii.crop,-pi/2+angles.crop,positions.crop_x0,-positions.crop_y0);
    
    tools.dbg('interval\t[%.3f\t%.3f],\tlength = %.3f,\tcrop2\n',interval_list(interval_ii).q(1),interval_list(interval_ii).q(2),arc_lengths.crop);
    interval_ii = interval_ii + 1;
    
    
end