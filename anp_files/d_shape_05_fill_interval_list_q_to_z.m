function interval_list = d_shape_05_fill_interval_list_q_to_z(in_params,in_data)

    if any(strcmp(who('global'),'debug_graphics'))
        global debug_graphics;
    else
        debug_graphics = false;
    end
    
    radii =                     in_params.radii;
    angles =                    in_params.angles;
    halfsecant_pole =           in_params.actual_halfsecant_pole;
    halfsecant_zero =           in_params.actual_halfsecant_zero;    
    im_pz_sorted =              in_data.im_pz_sorted;
    pole_zero_combinations =    in_data.im_pz_combinations;
    interval_list =             in_data.interval_list;
    
    secant_pole = 2*halfsecant_pole;
    secant_zero = 2*halfsecant_zero;
    [radii,arc_lengths] = calculate_detour_params(radii,angles,secant_pole,secant_zero);
    positions.crop_y0 = sqrt(radii.inf^2-2*radii.crop*radii.inf);
    positions.crop_x0 = radii.crop;
    angles.detour_pole_phi0 = asin(secant_pole/(2*radii.detour_pole));
    angles.detour_zero_phi0 = asin(secant_zero/(2*radii.detour_zero));
    
    % do the first element
    % can be
    % either: a negative imaginary pz with overlap up to a positive one that
    %         just touches the origin
    % or: the first straight interval
    if ~isempty(im_pz_sorted)
        im_pz_sorted = identify_border_cases(im_pz_sorted,halfsecant_pole,halfsecant_zero);
    end
    
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
    
    % ---------------------------------------------------------------------------------------------------------------------------------------------
    % treat positive p/z

    while(positive_pz_remain)
        idx_current_pz = idx_current_pz + 1;
        % there are positive poles/zeros left to treat
        % make a linear region to the next one
        interval_list(interval_ii).type = 'axis';
        interval_list(interval_ii).q(1) = prev_upper_bound;
        
        if isnan(idx_current_pz)
            idx_current_pz = idx_first_positive;
            interval_length = im_pz_sorted(idx_current_pz).value - halfsecant_pole*im_pz_sorted(idx_current_pz).pole - halfsecant_zero*im_pz_sorted(idx_current_pz).zero;
            za = 0;
        else
            interval_length = pole_zero_combinations(idx_current_pz-1).distance - halfsecant_pole*(im_pz_sorted(idx_current_pz-1).pole + im_pz_sorted(idx_current_pz).pole) - halfsecant_zero*(im_pz_sorted(idx_current_pz-1).zero + im_pz_sorted(idx_current_pz).zero);            
            za = im_pz_sorted(idx_current_pz-1).value + im_pz_sorted(idx_current_pz-1).pole*halfsecant_pole + im_pz_sorted(idx_current_pz-1).zero*halfsecant_zero;
        end
        interval_list(interval_ii).q_len = interval_length;
        
        interval_list(interval_ii).q(2) = interval_list(interval_ii).q(1) + interval_length;
        prev_upper_bound = interval_list(interval_ii).q(2);
        
        
        zb = im_pz_sorted(idx_current_pz).value - im_pz_sorted(idx_current_pz).pole*halfsecant_pole - im_pz_sorted(idx_current_pz).zero*halfsecant_zero;
        interval_list(interval_ii).input_fct_handle = @(q) im_axis_line(q,interval_list(interval_ii).q(1),interval_list(interval_ii).q(2),za,zb);
        
        tools.dbg('interval\t[%.3f\t%.3f],\tlength = %.3f,\tlinear_pos\n',interval_list(interval_ii).q(1),interval_list(interval_ii).q(2),interval_length);
        interval_ii = interval_ii + 1;
        
        % treat the pz
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
        
        tools.dbg('interval\t[%.3f\t%.3f],\tlength = %.3f,\tdetour_pos\n',interval_list(interval_ii).q(1),interval_list(interval_ii).q(2),interval_length);
        interval_ii = interval_ii + 1;
        
        positive_pz_remain = positive_pz_remain - 1;
    end
    
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
    
    % see p.41
    % 0
    if ~any([[im_pz_sorted.value] < 0,[im_pz_sorted.neg_overlapping],[im_pz_sorted.pos_overlapping],[im_pz_sorted.neg_on_origin]])
        % no neg pz, no pos overlap, already treated after crop2
    
    
    % 1
    elseif isnan(idx_current_pz)
            % have to add one last detour
            idx_current_pz = min([idx_last_negative,idx_first_positive]);
            interval_list(interval_ii).q(1) = prev_upper_bound;

            if im_pz_sorted(idx_current_pz).neg_on_origin
                interval_list(interval_ii).type = get_detour_type(im_pz_sorted(idx_current_pz).type,[]);
                interval_length = im_pz_sorted(idx_current_pz).pole*arc_lengths.pole + im_pz_sorted(idx_current_pz).zero*arc_lengths.pole;
                interval_list(interval_ii).q_len = interval_length;
                interval_list(interval_ii).q(2) = interval_list(interval_ii).q(1) + interval_length;

            elseif im_pz_sorted(idx_current_pz).neg_overlapping || im_pz_sorted(idx_current_pz).pos_overlapping
                interval_list(interval_ii).type = get_detour_type(im_pz_sorted(idx_current_pz).type,'part');
                arc_length_overlapping_pole = radii.detour_pole * (asin(-im_pz_sorted(idx_current_pz).value/radii.detour_pole) + angles.detour_pole_phi0);
                arc_length_overlapping_zero = radii.detour_zero * (asin(-im_pz_sorted(idx_current_pz).value/radii.detour_zero) + angles.detour_zero_phi0);
                interval_length = im_pz_sorted(idx_current_pz).pole*arc_length_overlapping_pole + im_pz_sorted(idx_current_pz).zero*arc_length_overlapping_zero;
                interval_list(interval_ii).q_len = interval_length;
                interval_list(interval_ii).q(2) = interval_list(interval_ii).q(1) + interval_length;

            else
                error('Oops, we shouldn''t be here. Apologies! Please report this crash to ricklis@student.ethz.ch together with the input you used.');
            end

            switch im_pz_sorted(idx_current_pz).type
                case 'p'
                    interval_list(interval_ii).input_fct_handle = @(q) circ_detour(map(q,interval_list(interval_ii).q(1),interval_list(interval_ii).q(2),0,interval_list(interval_ii).q_len),radii.detour_pole,secant_pole,im_pz_sorted(idx_current_pz).value);
                case 'z'
                    interval_list(interval_ii).input_fct_handle = @(q) circ_detour(map(q,interval_list(interval_ii).q(1),interval_list(interval_ii).q(2),0,interval_list(interval_ii).q_len),radii.detour_zero,secant_zero,im_pz_sorted(idx_current_pz).value);
            end
        
        tools.dbg('interval\t[%.3f\t%.3f],\tlength = %.3f,\tdetour\n',interval_list(interval_ii).q(1),interval_list(interval_ii).q(2),interval_length);
        
    % 2
    elseif ~isnan(idx_current_pz) && ~any([[im_pz_sorted.neg_overlapping],[im_pz_sorted.pos_overlapping],[im_pz_sorted.neg_on_origin]])
        % have to add one last linear interval to origin
            interval_list(interval_ii).type = 'axis';
            interval_list(interval_ii).q(1) = prev_upper_bound;
            
            interval_length = abs(im_pz_sorted(idx_current_pz).value + halfsecant_pole*im_pz_sorted(idx_current_pz).pole + halfsecant_zero*im_pz_sorted(idx_current_pz).zero);
            interval_list(interval_ii).q_len = interval_length;
            
            interval_list(interval_ii).q(2) = interval_list(interval_ii).q(1) + interval_length;
            
            za = im_pz_sorted(idx_current_pz).value + im_pz_sorted(idx_current_pz).pole*halfsecant_pole + im_pz_sorted(idx_current_pz).zero*halfsecant_zero;
            zb = 0;
            interval_list(interval_ii).input_fct_handle = @(q) im_axis_line(q,interval_list(interval_ii).q(1),interval_list(interval_ii).q(2),za,zb);
            
        tools.dbg('interval\t[%.3f\t%.3f],\tlength = %.3f,\tlinear_neg\n',interval_list(interval_ii).q(1),interval_list(interval_ii).q(2),interval_length);
        
        
    % 3
    elseif ~isnan(idx_current_pz) && any([[im_pz_sorted.neg_overlapping],[im_pz_sorted.pos_overlapping],[im_pz_sorted.neg_on_origin]])
        tools.dbg('negative and overlapping p/z\n');
            % lin interval first
            interval_list(interval_ii).type = 'axis';
            interval_list(interval_ii).q(1) = prev_upper_bound;
            
            interval_length = pole_zero_combinations(idx_current_pz).distance - halfsecant_pole*(im_pz_sorted(idx_current_pz).pole + im_pz_sorted(idx_current_pz+1).pole) - halfsecant_zero*(im_pz_sorted(idx_current_pz).zero + im_pz_sorted(idx_current_pz+1).zero);            
            interval_list(interval_ii).q_len = interval_length;
            
            interval_list(interval_ii).q(2) = interval_list(interval_ii).q(1) + interval_length;
            prev_upper_bound = interval_list(interval_ii).q(2);
            
            za = im_pz_sorted(idx_current_pz).value   + im_pz_sorted(idx_current_pz).pole  *halfsecant_pole + im_pz_sorted(idx_current_pz).zero  *halfsecant_zero;
            zb = im_pz_sorted(idx_current_pz+1).value - im_pz_sorted(idx_current_pz+1).pole*halfsecant_pole - im_pz_sorted(idx_current_pz+1).zero*halfsecant_zero;
            interval_list(interval_ii).input_fct_handle = @(q) im_axis_line(q,interval_list(interval_ii).q(1),interval_list(interval_ii).q(2),za,zb);

        tools.dbg('interval\t[%.3f\t%.3f],\tlength = %.3f,\tlinear_neg\n',interval_list(interval_ii).q(1),interval_list(interval_ii).q(2),interval_length);
            
            % do detour here
            interval_ii = interval_ii + 1;
            idx_current_pz = idx_current_pz + 1;
            interval_list(interval_ii).type = get_detour_type(im_pz_sorted(idx_current_pz).type,'part');
            interval_list(interval_ii).q(1) = prev_upper_bound;
            
            arc_length_overlapping_pole = radii.detour_pole * (asin(-im_pz_sorted(idx_current_pz).value/radii.detour_pole) + angles.detour_pole_phi0);
            arc_length_overlapping_zero = radii.detour_zero * (asin(-im_pz_sorted(idx_current_pz).value/radii.detour_zero) + angles.detour_zero_phi0);
            interval_length = im_pz_sorted(idx_current_pz).pole*(arc_lengths.detour_pole - arc_length_overlapping_pole) + im_pz_sorted(idx_current_pz).zero*(arc_lengths.detour_zero - arc_length_overlapping_zero);
            interval_list(interval_ii).q_len = interval_length;
            
            interval_list(interval_ii).q(2) = interval_list(interval_ii).q(1) + interval_length;
            
            switch im_pz_sorted(idx_current_pz).type
                case 'p'
                    interval_list(interval_ii).input_fct_handle = @(q) circ_detour(map(q,interval_list(interval_ii).q(1),interval_list(interval_ii).q(2),0,interval_list(interval_ii).q_len),radii.detour_pole,secant_pole,im_pz_sorted(idx_current_pz).value);
                case 'z'
                    interval_list(interval_ii).input_fct_handle = @(q) circ_detour(map(q,interval_list(interval_ii).q(1),interval_list(interval_ii).q(2),0,interval_list(interval_ii).q_len),radii.detour_zero,secant_zero,im_pz_sorted(idx_current_pz).value);
            end
            
            
        tools.dbg('interval\t[%.3f\t%.3f],\tlength = %.3f,\tdetour\n',interval_list(1).q(1),interval_list(1).q(2),interval_length);
    else
        error('Oops, we shouldn''t be here. Apologies! Please report this crash to ricklis@student.ethz.ch together with the input you used.');
    end
    
    
    if debug_graphics
        figure;
        for kk = 1:length(interval_list)
            t = intvl(interval_list(kk).q);
            scatter(real(interval_list(kk).input_fct_handle(t)),imag(interval_list(kk).input_fct_handle(t)));
            hold on;
        end
        axis equal;
    end
    
end

function pz_list_sorted = identify_border_cases(pz_list_sorted,halfsecant_pole,halfsecant_zero)
    for ii = 1:length(pz_list_sorted)
        pz_value = pz_list_sorted(ii).value;
        switch(pz_list_sorted(ii).type)
            case 'p'
                if pz_value < 0 && (pz_value + halfsecant_pole == 0)
                    pz_list_sorted(ii).neg_on_origin = true;
                elseif pz_value < 0 && (pz_value + halfsecant_pole > 0)
                    pz_list_sorted(ii).neg_overlapping = true;
                elseif pz_value >= 0 && (pz_value - halfsecant_pole == 0)
                    pz_list_sorted(ii).pos_on_origin = true;
                elseif pz_value >= 0 && (pz_value - halfsecant_pole < 0)
                    pz_list_sorted(ii).pos_overlapping = true;
                end
            case 'z'
                if pz_value < 0 && (pz_value + halfsecant_zero == 0)
                    pz_list_sorted(ii).neg_on_origin = true;
                elseif pz_value < 0 && (pz_value + halfsecant_zero > 0)
                    pz_list_sorted(ii).neg_overlapping = true;
                elseif pz_value >= 0 && (pz_value - halfsecant_zero == 0)
                    pz_list_sorted(ii).pos_on_origin = true;
                elseif pz_value >= 0 && (pz_value - halfsecant_zero < 0)
                    pz_list_sorted(ii).pos_overlapping = true;
                end
            otherwise
                tools.dbg('identify_border_cases: ii = %d, type = %s, value = %f\n',ii,pz_list_sorted(ii).type,pz_list_sorted(ii).value);
                error('Oops, we shouldn''t be here. Apologies! Please report this crash to ricklis@student.ethz.ch together with the input you used.');
        end
    end
end

function type = get_detour_type(pz_type,remark)
    underline = [];
    if ~isempty(remark)
        underline = '_';
    end
    
    switch pz_type
        case 'p'
            type = ['detour_pole',underline,remark];
        case 'z'
            type = ['detour_zero',underline,remark];
        otherwise
            error('Oops, we shouldn''t be here. Apologies! Please report this crash to ricklis@student.ethz.ch together with the input you used.');
    end 
end

function [radii,arc_lengths] = calculate_detour_params(radii,angles,secant_pole,secant_zero)
    radii.crop = radii.inf*sin(angles.crop)/(1+sin(angles.crop));
    
    radii.detour_pole = secant_pole*sqrt(4+tan(angles.detour)^2)/4;
    radii.detour_zero = secant_zero*sqrt(4+tan(angles.detour)^2)/4;
    
    arc_lengths.inf = (pi-2*angles.crop)*radii.inf;
    arc_lengths.crop = (pi/2+angles.crop)*radii.crop;
    arc_lengths.detour_pole = 2*radii.detour_pole*asin(secant_pole/(2*radii.detour_pole));
    arc_lengths.detour_zero = 2*radii.detour_zero*asin(secant_zero/(2*radii.detour_zero));
end

function z = circ_normal(q,R,phi_0,x0,y0)
    assert((real(y0) == 0 && imag(y0) ~= 0) || (real(y0) ~= 0 && imag(y0) == 0) || (y0 == 0));
    y0 = real(y0) + imag(y0);
    
    z = R*exp(1i*(phi_0 + q/R)) + x0 + 1i*y0;
end

function z = circ_detour(q,R,sec,y_pz)
    y_pz = real(y_pz) + imag(y_pz);
    
    z = R*(exp(1i*(q/R - asin(sec/(2*R)))) - sqrt(1-(sec/(2*R))^2)) + 1i*y_pz;
end

function z = im_axis_line(q,qa,qb,za,zb)
    z = 1i*map(q,qa,qb,za,zb);
end

function y = map(x,t0,t1,u0,u1)
    y = ((u0-u1).*x + (t0*u1-t1*u0))/(t0-t1);
end

function t = intvl(in)
    t = linspace(in(1),in(2),100);
end