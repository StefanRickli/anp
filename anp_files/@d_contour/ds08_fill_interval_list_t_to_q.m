function [] = ds08_fill_interval_list_t_to_q(this)
    
    interval_list =     this.interval_list;
    
    this.ds08a_calc_shares_of_t_intervals();
    
    n_axis_intvls =     sum(strncmp({interval_list.type},'axis',4));
    n_pole_detours =    sum(strncmp({interval_list.type},'detour_pole',11)) - sum(strncmp({interval_list.type},'detour_pole_part',16))/2;
    n_zero_detours =    sum(strncmp({interval_list.type},'detour_zero',11)) - sum(strncmp({interval_list.type},'detour_zero_part',16))/2;
    
    for ii = 1:length(interval_list)
        switch interval_list(ii).type
            case 'axis'
                interval_list(ii).t_len = this.shares.axis/n_axis_intvls;
            case 'crop'
                interval_list(ii).t_len = this.shares.crop;
            case 'inf'
                interval_list(ii).t_len = this.shares.inf;
            case 'detour_pole_part'
                interval_list(ii).t_len = this.shares.poles/(2*n_pole_detours);
            case 'detour_pole'
                interval_list(ii).t_len = this.shares.poles/n_pole_detours;
            case 'detour_zero_part'
                interval_list(ii).t_len = this.shares.zeros/(2*n_zero_detours);
            case 'detour_zero'
                interval_list(ii).t_len = this.shares.zeros/n_zero_detours;
            otherwise
                error('Oops, we shouldn''t be here. Apologies! Please report this crash to stefanrickli [at] gmx.ch together with the input you used.');
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
                ii_prev = one_indexing_modulo(ii-1,length(interval_list));
                
                ta_prev = interval_list(ii_prev).t(1);
                tb_prev = interval_list(ii_prev).t(2);
                qa_prev = interval_list(ii_prev).q(1);
                qb_prev = interval_list(ii_prev).q(2);
                va = (qb_prev - qa_prev)/(tb_prev - ta_prev);
            end
            
            if ii == length(interval_list)
                vb = 0.5;
            else
                ii_next = one_indexing_modulo(ii + 1,length(interval_list));

                ta_next = interval_list(ii_next).t(1);
                tb_next = interval_list(ii_next).t(2);
                qa_next = interval_list(ii_next).q(1);
                qb_next = interval_list(ii_next).q(2);
                vb = (qb_next - qa_next)/(tb_next - ta_next);
            end
            
            % Calculate the parameters for the nonlinear mapping function
            % 'exp_Lorentz' for this interval. We give it the starting
            % point (ta,qa), end point (tb,qb) and a desired slope at the
            % beginning va and at the end vb.
            [c_lorentz,gamma_star,r1,r2,c1,c2] =    this.ds09_mixed_exp_lorentz(ta,tb,qa,qb,va,vb);
            
            % The solver in the interpolation function delivers not so
            % precise start and endpoints. We need to align them perfectly
            % with the desired qa and qb, so we simply do a linear mapping.
            qa_real =                               exp_Lorentz(ta,ta,tb,qa,c_lorentz,gamma_star,r1,r2,c1,c2);
            qb_real =                               exp_Lorentz(tb,ta,tb,qa,c_lorentz,gamma_star,r1,r2,c1,c2);
            interval_list(ii).density_fct_handle =  @(t) map(exp_Lorentz(t,ta,tb,qa,c_lorentz,gamma_star,r1,r2,c1,c2), qa_real,qb_real,qa,qb);
            
        elseif any(strcmp(interval_list(ii).type,{'detour_pole','detour_zero','detour_pole_part','detour_zero_part','crop','inf'}))
            interval_list(ii).density_fct_handle =  @(t) map(t,ta,tb,qa,qb);
        else
            error('Oops, we shouldn''t be here. Apologies! Please report this crash to stefanrickli [at] gmx.ch together with the input you used.');
        end
    end
    
    this.interval_list = interval_list;
end
