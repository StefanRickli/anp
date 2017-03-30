classdef d_contour < handle
    properties(SetAccess = private)
        t_values
        z_values
        
        poles
        zeros
        im_pz_sorted
        im_pz_combinations
        interval_list
        
        radii
        angles
        secant_pole
        secant_zero
        halfsecant_pole
        halfsecant_pole_max
        halfsecant_zero
        halfsecant_zero_max
        halfsecant_margin
        arc_lengths
        positions
        weight
        shares
    end
    
    methods(Access = public)
        function this = d_contour(t,in_params)
            this.t_values = t;
            
            this.poles =    in_params.poles;
            this.zeros =    in_params.zeros;
            this.radii =    in_params.radii;
            this.angles =   in_params.angles;
            this.halfsecant_pole_max =  in_params.halfsecant_pole_max;
            this.halfsecant_zero_max =  in_params.halfsecant_zero_max;
            this.halfsecant_margin =    in_params.halfsecant_margin;
            this.weight =   in_params.weight;
            
            this.run();
        end
        
        function res = get_z_values(this)
            res = this.z_values;
        end
    end
    
    methods(Access = private)
        function [] = run(this)
            this.ds01_sort_poles_and_zeros();
            this.ds02_list_pole_and_zero_combinations();
            this.ds03_calc_definitive_halfsecants();
            this.ds04_init_interval_list();
            this.ds05_calc_radii_arcs_angles_positions();
            this.ds06_identify_border_cases();
            this.ds07_fill_interval_list_q_to_z();
            this.ds08_fill_interval_list_t_to_q();
            this.evaluate();
        end
                
        function ds07_fill_interval_list_q_to_z(this)
            if any(strcmp(who('global'),'debug_graphics'))
                global debug_graphics;
            else
                debug_graphics = false;
            end

            [interval_ii,idx_current_pz,prev_upper_bound,positive_pz_remain,idx_first_positive] =  this.ds07a_1st_interval();
            [interval_ii,idx_current_pz,prev_upper_bound] =                     this.ds07b_positive_pz(interval_ii,idx_current_pz,prev_upper_bound,positive_pz_remain,idx_first_positive);
            [interval_ii] =                                                     this.ds07c_last_straight_before_crop1(interval_ii,idx_current_pz,prev_upper_bound);

            [interval_ii] = this.ds07d_crop1_inf_crop2(interval_ii);

            [interval_ii,idx_current_pz,prev_upper_bound] =                     this.ds07g_first_after_crop2(interval_ii);
            [interval_ii,idx_current_pz,prev_upper_bound] =                     this.ds07h_negative_pz(interval_ii,idx_current_pz,prev_upper_bound);
                                                                                this.ds07i_last_interval(interval_ii,idx_current_pz,prev_upper_bound);

            if debug_graphics
                figure;
                for kk = 1:length(this.interval_list)
                    t = intvl(this.interval_list(kk).q);
                    scatter(real(this.interval_list(kk).input_fct_handle(t)),imag(this.interval_list(kk).input_fct_handle(t)));
                    hold on;
                end
                axis equal;
            end
        end
        
        function [] = evaluate(this)
            % Performs the mapping t |--> q |--> z for every point in the vector t
            % 
            % PRE: Monotonically increasing list of real values in the range [0,1]
            
            % As we know that the values in t are monotonically increasing, there's
            % no need to look for the right interval for every point of t. Simply
            % check whether t(ii) > current interval's upper limit and switch to
            % the next interval if the comparison yields true.
            ii_interval = 1;
            this.z_values = zeros(size(this.t_values));
            for ii = 1:length(this.t_values)
                while (ii_interval + 1 <= length(this.interval_list)) && (this.t_values(ii) > this.interval_list(ii_interval).t(2))
                    ii_interval = ii_interval + 1;
                end
                
                % Map t(ii) |--> q(ii)
                Q = this.interval_list(ii_interval).density_fct_handle(this.t_values(ii));
                
                % Map q(ii) |--> z(ii)
                this.z_values(ii) = this.interval_list(ii_interval).input_fct_handle(Q);
            end
        end
        
    end
end