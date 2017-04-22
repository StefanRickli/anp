function [] = ds07g_last_interval(this,interval_ii,idx_current_pz,prev_upper_bound,idx_first_positive,idx_last_negative)
    im_pz_sorted =          this.im_pz_sorted;
    interval_list =         this.interval_list;
    secant_pole =           this.secant_pole;
    secant_zero =           this.secant_zero;
    halfsecant_pole =       this.halfsecant_pole;
    halfsecant_zero =       this.halfsecant_zero;
    radii =                 this.radii;
    pole_zero_combinations =this.im_pz_combinations;
    angles =                this.angles;
    arc_lengths =           this.arc_lengths;
    
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
                interval_list(interval_ii).type = [repmat('detour_pole', im_pz_sorted(idx_current_pz).pole), ...
                                                   repmat('detour_zero', im_pz_sorted(idx_current_pz).zero)];
                interval_length = arc_lengths.detour_pole * im_pz_sorted(idx_current_pz).pole + ...
                                  arc_lengths.detour_zero * im_pz_sorted(idx_current_pz).zero;
                interval_list(interval_ii).q_len = interval_length;
                interval_list(interval_ii).q(2) = interval_list(interval_ii).q(1) + interval_length;

            elseif im_pz_sorted(idx_current_pz).neg_overlapping || im_pz_sorted(idx_current_pz).pos_overlapping
                interval_list(interval_ii).type = [repmat('detour_pole_part', im_pz_sorted(idx_current_pz).pole), ...
                                                   repmat('detour_zero_part', im_pz_sorted(idx_current_pz).zero)];
                arc_length_overlapping_pole = radii.detour_pole * ( asin(-im_pz_sorted(idx_current_pz).value/radii.detour_pole) + ...
                                                                    angles.detour_pole_phi0 ...
                                                                   );
                arc_length_overlapping_zero = radii.detour_zero * ( asin(-im_pz_sorted(idx_current_pz).value/radii.detour_zero) + ...
                                                                    angles.detour_zero_phi0 ...
                                                                   );
                interval_length = arc_length_overlapping_pole * im_pz_sorted(idx_current_pz).pole + ...
                                  arc_length_overlapping_zero * im_pz_sorted(idx_current_pz).zero;
                interval_list(interval_ii).q_len = interval_length;
                interval_list(interval_ii).q(2) = interval_list(interval_ii).q(1) + interval_length;

            else
                error('Oops, we shouldn''t be here. Apologies! Please report this crash to stefanrickli [at] gmx.ch together with the input you used.');
            end

            switch im_pz_sorted(idx_current_pz).type
                case 'p'
                    interval_list(interval_ii).input_fct_handle = @(q) circ_detour(map(q,interval_list(interval_ii).q(1),interval_list(interval_ii).q(2),0,interval_list(interval_ii).q_len),radii.detour_pole,secant_pole,im_pz_sorted(idx_current_pz).value);
                case 'z'
                    interval_list(interval_ii).input_fct_handle = @(q) circ_detour(map(q,interval_list(interval_ii).q(1),interval_list(interval_ii).q(2),0,interval_list(interval_ii).q_len),radii.detour_zero,secant_zero,im_pz_sorted(idx_current_pz).value);
            end
        
        dbg_out('interval\t[%.3f\t%.3f],\tlength = %.3f,\tdetour\n',interval_list(interval_ii).q(1),interval_list(interval_ii).q(2),interval_length);
        
    % 2
    elseif ~isnan(idx_current_pz) && ~any([[im_pz_sorted.neg_overlapping],[im_pz_sorted.pos_overlapping],[im_pz_sorted.neg_on_origin]])
        % have to add one last linear interval to origin
            interval_list(interval_ii).type = 'axis';
            interval_list(interval_ii).q(1) = prev_upper_bound;
            
            interval_length = abs( im_pz_sorted(idx_current_pz).value + ...
                                   halfsecant_pole * im_pz_sorted(idx_current_pz).pole + ...
                                   halfsecant_zero * im_pz_sorted(idx_current_pz).zero ...
                                  );
            interval_list(interval_ii).q_len = interval_length;
            
            interval_list(interval_ii).q(2) = interval_list(interval_ii).q(1) + interval_length;
            
            za = im_pz_sorted(idx_current_pz).value + ...
                 halfsecant_pole * im_pz_sorted(idx_current_pz).pole + ...
                 halfsecant_zero * im_pz_sorted(idx_current_pz).zero;
            zb = 0;
            interval_list(interval_ii).input_fct_handle = @(q) im_axis_line(q,interval_list(interval_ii).q(1),interval_list(interval_ii).q(2),za,zb);
            
        dbg_out('interval\t[%.3f\t%.3f],\tlength = %.3f,\tlinear_neg\n',interval_list(interval_ii).q(1),interval_list(interval_ii).q(2),interval_length);
        
        
    % 3
    elseif ~isnan(idx_current_pz) && any([[im_pz_sorted.neg_overlapping],[im_pz_sorted.pos_overlapping],[im_pz_sorted.neg_on_origin]])
        dbg_out('negative and overlapping p/z\n');
            % lin interval first
            interval_list(interval_ii).type = 'axis';
            interval_list(interval_ii).q(1) = prev_upper_bound;
            
            interval_length = pole_zero_combinations(idx_current_pz).distance + ...
                              -halfsecant_pole * ( im_pz_sorted(idx_current_pz).pole + ...
                                                   im_pz_sorted(idx_current_pz+1).pole) + ...
                              -halfsecant_zero * ( im_pz_sorted(idx_current_pz).zero + ...
                                                   im_pz_sorted(idx_current_pz+1).zero);
            interval_list(interval_ii).q_len = interval_length;
            
            interval_list(interval_ii).q(2) = interval_list(interval_ii).q(1) + interval_length;
            prev_upper_bound = interval_list(interval_ii).q(2);
            
            za = im_pz_sorted(idx_current_pz).value + ...
                 halfsecant_pole * im_pz_sorted(idx_current_pz).pole + ...
                 halfsecant_zero * im_pz_sorted(idx_current_pz).zero;
            zb = im_pz_sorted(idx_current_pz+1).value + ...
                 -halfsecant_pole * im_pz_sorted(idx_current_pz+1).pole + ...
                 -halfsecant_zero * im_pz_sorted(idx_current_pz+1).zero;
            interval_list(interval_ii).input_fct_handle = @(q) im_axis_line(q,interval_list(interval_ii).q(1),interval_list(interval_ii).q(2),za,zb);

        dbg_out('interval\t[%.3f\t%.3f],\tlength = %.3f,\tlinear_neg\n',interval_list(interval_ii).q(1),interval_list(interval_ii).q(2),interval_length);
            
            % do detour here
            interval_ii = interval_ii + 1;
            idx_current_pz = idx_current_pz + 1;
            interval_list(interval_ii).type = [repmat('detour_pole_part', im_pz_sorted(idx_current_pz).pole), ...
                                               repmat('detour_zero_part', im_pz_sorted(idx_current_pz).zero)];
            interval_list(interval_ii).q(1) = prev_upper_bound;
            
            arc_length_overlapping_pole = radii.detour_pole * ( asin(-im_pz_sorted(idx_current_pz).value/radii.detour_pole) + ...
                                                                angles.detour_pole_phi0 ...
                                                               );
            arc_length_overlapping_zero = radii.detour_zero * ( asin(-im_pz_sorted(idx_current_pz).value/radii.detour_zero) + ...
                                                                angles.detour_zero_phi0 ...
                                                               );
            interval_length = (arc_lengths.detour_pole - arc_length_overlapping_pole) * im_pz_sorted(idx_current_pz).pole + ...
                              (arc_lengths.detour_zero - arc_length_overlapping_zero) * im_pz_sorted(idx_current_pz).zero;
            interval_list(interval_ii).q_len = interval_length;
            
            interval_list(interval_ii).q(2) = interval_list(interval_ii).q(1) + interval_length;
            
            switch im_pz_sorted(idx_current_pz).type
                case 'p'
                    interval_list(interval_ii).input_fct_handle = @(q) circ_detour(map(q,interval_list(interval_ii).q(1),interval_list(interval_ii).q(2),0,interval_list(interval_ii).q_len), ...
                                                                                   radii.detour_pole, ...
                                                                                   secant_pole, ...
                                                                                   im_pz_sorted(idx_current_pz).value);
                case 'z'
                    interval_list(interval_ii).input_fct_handle = @(q) circ_detour(map(q,interval_list(interval_ii).q(1),interval_list(interval_ii).q(2),0,interval_list(interval_ii).q_len), ...
                                                                                   radii.detour_zero, ...
                                                                                   secant_zero, ...
                                                                                   im_pz_sorted(idx_current_pz).value);
            end
            
            
        dbg_out('interval\t[%.3f\t%.3f],\tlength = %.3f,\tdetour\n',interval_list(1).q(1),interval_list(1).q(2),interval_length);
    else
        error('Oops, we shouldn''t be here. Apologies! Please report this crash to stefanrickli [at] gmx.ch together with the input you used.');
    end
    
    this.interval_list = interval_list;    
end