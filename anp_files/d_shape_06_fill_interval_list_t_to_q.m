function interval_list = d_shape_06_fill_interval_list_t_to_q(in_params,in_data)
    
    poles =         in_params.poles;
    zeros =         in_params.zeros;
    im_pz_sorted =  in_data.im_pz_sorted;
    interval_list = in_data.interval_list;
    
    n_im_poles_distinct = sum([im_pz_sorted.type] == 'p');
    n_im_zeros_distinct = sum([im_pz_sorted.type] == 'z');
    tf_relative_degree = length(poles) - length(zeros);
    shares = get_shares(tf_relative_degree,n_im_poles_distinct,n_im_zeros_distinct);
    shares = determine_crop_inf_shares(shares,interval_list);
 
    n_axis_intvls = sum(strncmp({interval_list.type},'axis',4));
    n_pole_detours = sum(strncmp({interval_list.type},'detour_pole',11)) - sum(strncmp({interval_list.type},'detour_pole_part',16))/2;
    n_zero_detours = sum(strncmp({interval_list.type},'detour_zero',11)) - sum(strncmp({interval_list.type},'detour_zero_part',16))/2;
    
    for ii = 1:length(interval_list)
        switch interval_list(ii).type
            case 'axis'
                interval_list(ii).t_len = shares.axis/n_axis_intvls;
            case 'crop'
                interval_list(ii).t_len = shares.crop;
            case 'inf'
                interval_list(ii).t_len = shares.inf;
            case 'detour_pole_part'
                interval_list(ii).t_len = shares.poles/(2*n_pole_detours);
            case 'detour_pole'
                interval_list(ii).t_len = shares.poles/n_pole_detours;
            case 'detour_zero_part'
                interval_list(ii).t_len = shares.zeros/(2*n_zero_detours);
            case 'detour_zero'
                interval_list(ii).t_len = shares.zeros/n_zero_detours;
            otherwise
                error('Oops, we shouldn''t be here. Apologies! Please report this crash to ricklis@student.ethz.ch together with the input you used.');
        end
    end
    
    interval_list(1).t(1) = 0;
    interval_list(1).t(2) = interval_list(1).t_len;
    for ii = 2:length(interval_list)
        interval_list(ii).t(1) = interval_list(ii-1).t(2);
        interval_list(ii).t(2) = interval_list(ii).t(1) + interval_list(ii).t_len;
    end
    
    for ii = 1:length(interval_list)
        ta = interval_list(ii).t(1);
        tb = interval_list(ii).t(2);
        qa = interval_list(ii).q(1);
        qb = interval_list(ii).q(2);

        if strcmp(interval_list(ii).type,'axis')            
            if ii == 1
                va = 0.5;
            else
                ii_prev = iterator_modulo(ii-1,length(interval_list));
                
                ta_prev = interval_list(ii_prev).t(1);
                tb_prev = interval_list(ii_prev).t(2);
                qa_prev = interval_list(ii_prev).q(1);
                qb_prev = interval_list(ii_prev).q(2);
                va = (qb_prev - qa_prev)/(tb_prev - ta_prev);
            end
            
            if ii == length(interval_list)
                vb = 0.5;
            else
                ii_next = iterator_modulo(ii + 1,length(interval_list));

                ta_next = interval_list(ii_next).t(1);
                tb_next = interval_list(ii_next).t(2);
                qa_next = interval_list(ii_next).q(1);
                qb_next = interval_list(ii_next).q(2);
                vb = (qb_next - qa_next)/(tb_next - ta_next);
            end
            
            [c_lorentz,gamma_star,r1,r2,c1,c2] = d_shape_07_mixed_exp_lorentz(ta,tb,qa,qb,va,vb);
            qa_real = exp_Lorentz(ta,ta,tb,qa,c_lorentz,gamma_star,r1,r2,c1,c2);
            qb_real = exp_Lorentz(tb,ta,tb,qa,c_lorentz,gamma_star,r1,r2,c1,c2);
            interval_list(ii).density_fct_handle = @(t) map(exp_Lorentz(t,ta,tb,qa,c_lorentz,gamma_star,r1,r2,c1,c2),qa_real,qb_real,qa,qb);
            
        elseif any(strcmp(interval_list(ii).type,{'detour_pole','detour_zero','detour_pole_part','detour_zero_part','crop','inf'}))
            interval_list(ii).density_fct_handle = @(t) map(t,ta,tb,qa,qb);
        else
            error('Oops, we shouldn''t be here. Apologies! Please report this crash to ricklis@student.ethz.ch together with the input you used.');
        end
    end
end

function shares = get_shares(tf_relative_degree,n_im_poles_distinct,n_im_zeros_distinct)
    assert(tf_relative_degree >= 0);
    
    weights.crop_inf = 3/sqrt(tf_relative_degree+1);
    weights.interpolation = 9;
    weights.pole = n_im_poles_distinct*1;
    weights.zero = n_im_zeros_distinct*1/3;
    temp = struct2cell(weights);
    weight_sum = sum([temp{:}]);
    
    shares.crop_inf = weights.crop_inf / weight_sum;
    shares.axis = weights.interpolation / weight_sum;
    shares.poles = weights.pole / weight_sum;
    shares.zeros = weights.zero / weight_sum;    
end

function shares = determine_crop_inf_shares(shares,interval_list)
    crop_idx = strncmp({interval_list.type},'crop',4);
    crop_idx = find(crop_idx == true,1,'first');
    crop_arc_length = interval_list(crop_idx).q_len;
    
    inf_idx = strncmp({interval_list.type},'inf',3);
    inf_idx = find(inf_idx == true,1,'first');
    inf_arc_length = interval_list(inf_idx).q_len;
    
    crop_inf_arc_length = 2*crop_arc_length + inf_arc_length;
    
    shares.crop = crop_arc_length/crop_inf_arc_length * shares.crop_inf;
    shares.inf =  inf_arc_length/crop_inf_arc_length * shares.crop_inf;
end

function y = map(x,t0,t1,u0,u1)
    y = ((u0-u1).*x + (t0*u1-t1*u0))/(t0-t1);
end

function y = iterator_modulo(x,m)
    y = mod(x - 1,m) + 1;
end

function z = exp_Lorentz(T,a,b,ya,c_lorentz,gamma_star,r1,r2,c1,c2)
    delta_x = b - a;
    
    if ~isnan(c_lorentz)
        Lorentz_0 = c_lorentz  *  (atan( (2*T-a-b) ./ (2*gamma_star) ) + ...
                                   atan(   delta_x / (2*gamma_star) ) + ...
                                   (4*gamma_star.*(a-T)) ./ (delta_x^(2)+4*gamma_star^2));
    else
        Lorentz_0 = zeros(size(T));            
    end

    G1 = c1 * (a*exp(-r1*b) + exp(-r1*a)/r1 - exp(-r1*T)/r1 - T*exp(-r1*b));
    G2 = c2 * (a*exp( r2*a) - exp( r2*a)/r2 + exp( r2*T)/r2 - T*exp( r2*a));
    G = G1 + G2;

    z = ya+G+Lorentz_0;

end
