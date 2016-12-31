% as in p.41
classdef anp_tf_processor < handle
    properties(SetAccess = private)
        % g: general
        g_uid
        g_wait          logical
        g_busy          logical
        
        % tf: transfer function properties
        tf_obj
        tf_poles        double
        tf_zeros        double
        tf_G
        
        % p: parameters
        p_time
        p_radii
        p_angles
        p_weights
        
        % d: calculated data
        d_time_points   double
        d_data_points   double
        d_z_values      double
        d_w_values      double
    end
    methods(Access = public)
        function this = anp_tf_processor()
            this.g_uid =        java.util.UUID.randomUUID.toString;
            this.g_wait =       false;
            this.g_busy =       false;
            this.d_z_values =  [];
            this.d_w_values = [];
            
            load('demo_values.mat');
            this.d_data_points = t;
            this.d_z_values = in_values;
            this.d_w_values = out_values;
            
            fprintf('anp_tf_processor\t%s: Instance created.\n',this.g_uid);
        end
        
        function [] = set_all_params(this,new_params)
            this.g_wait =       true;
            if any([this.set_tf(new_params.tf_obj),...
                    this.set_radii(new_params.radii),...
                    this.set_angles(new_params.angles),...
                    this.set_weights(new_params.weights)])
                this.recalculate();
            end
            this.g_wait =       false;
        end
        
        function dirty = set_tf(this,new_tf)
            if ~isa(new_tf,'tf')
                error('Expected an argument of type ''tf''.');
            end
            if ~isequal(this.tf_obj,new_tf)
                this.tf_obj =   new_tf;
                this.tf_poles = roots(new_tf.Denominator{1})';
                this.tf_zeros = roots(new_tf.Numerator{1})';
                this.tf_G =     @(z) polyval(poly(zeros),z)./polyval(poly(poles),z);
                
                dirty =         true;
                if ~this.g_wait
                    this.recalculate();
                end
            else
                dirty =         false;
                disp('TF: Nothing changed.');
            end
            
            this.echo_tf();
        end
        
        function dirty = set_radii(this,new_radii)
            % TODO autocalc changes radii, so we can't really compare for
            % equality directly...
%             if ~(isvector(new_radii) && isnumeric(new_radii))
%                 error('Expected an array of complex values.');
%             end
            if isempty(this.p_radii) || ~all([this.p_radii.R,new_radii.R])
                this.p_radii =  new_radii;
                dirty =         true;
                if ~this.g_wait
                    this.recalculate();
                end
            else
                dirty =         false;
                disp('Radii: Nothing changed.');
            end
        end
        
        function dirty = set_angles(this,new_angles)
%             if ~(isvector(new_angles) && isnumeric(new_angles))
%                 error('Expected an array of complex values.');
%             end
            if isempty(this.p_angles) || ~all([this.p_angles.crop_inf_transition == new_angles.crop_inf_transition,...
                     this.p_angles.min_angle_contribution_at_R == new_angles.min_angle_contribution_at_R...
                     this.p_angles.detour == new_angles.detour])
                dirty =         true;
                this.p_angles =  new_angles;
                if ~this.g_wait
                    this.recalculate();
                end
            else
                dirty =         false;
                disp('Angles: Nothing changed.');
            end
        end
        
        function dirty = set_weights(this,new_weights)
%             if ~(isvector(new_weights) && isnumeric(new_weights))
%                 error('Expected an array of complex values.');
%             end
            if ~isequal(this.p_radii,new_weights)
                dirty =             true;
                this.p_weights =    new_weights;
                if ~this.g_wait
                    this.recalculate();
                end
            else
                dirty =         false;
                disp('Weights: Nothing changed.');
            end
        end
        
        % put the transfer function to the text output
        function [] = echo_tf(this)
            transfer_function = this.tf_obj %#ok<NOPRT,NASGU>
            fprintf('with\n');
            poles = this.tf_poles           %#ok<NOPRT,NASGU>
            fprintf('and\n');
            zeros = this.tf_zeros           %#ok<NOPRT,NASGU>

            if any([this.tf_obj.IODelay,this.tf_obj.InputDelay,this.tf_obj.OutputDelay])
                warning('Delay has been provided, but will be ignored, as this is currently unsupported.');
            end
        end
        
        function data = get_time_points(this)
            while this.g_busy % TODO potentially unsafe parallel processing
                pause(1/50);
            end
            
            data.time_points =   this.d_time_points;
            data.data_points =   this.d_data_points;
            data.time_props =    this.p_time;
        end
        
        function data = get_data(this)
            while this.g_busy % TODO potentially unsafe parallel processing
                pause(1/50);
            end
            
            % TODO assign more data fields according to what GUI expects
            
            data.z_values =     this.d_z_values;
            data.w_values =     this.d_w_values;
        end
        function delete(this)
            fprintf('anp_tf_processor\t%s: Deletion requested.\n',this.g_uid);
        end
    end
    methods(Access = private)
        function [] = recalculate(this)
            this.g_busy =       true;
            
            if this.p_radii.auto_main_R
                this.p_radii.R =    anp_calc_main_R(this.tf_poles,this.tf_zeros,this.p_angles.min_angle_contribution_at_R);
            end
            
            this.calc_time_params();
            
            this.g_busy =   false;
            notify(this,'tf_processor_data_ready');
        end
        
        function [] = calc_time_params(this)
            this.p_time.n_time_steps =          120; % TODO include this in argument parsing, don't hardcode!
            
            p = this.tf_poles;
            z = this.tf_zeros;
                p_neg_real =    sum(real(p) <= 0); % TODO check if equality needed
                p_pos_real =    sum(real(p) > 0);
                p_pure_imag =   sum(real(p) == 0);
                z_neg_real =    sum(real(z) <= 0); % TODO
                z_pos_real =    sum(real(z) > 0);
                z_pure_imag =   sum(real(z) == 0);
            
            total_tf_rel_deg =                  p_neg_real - z_neg_real - p_pos_real + z_pos_real;
            this.p_time.oversampling_factor =   ceil(total_tf_rel_deg + this.p_weights.pole * p_pure_imag + this.p_weights.zero * z_pure_imag);
            this.p_time.n_data_points =         this.p_time.n_time_steps * this.p_time.oversampling_factor;
            
            this.d_time_points =    linspace(0,1,this.p_time.n_time_steps+1);
            this.d_time_points =    this.d_time_points(1:end-1);
            this.d_data_points =    linspace(0,1,this.p_time.n_data_points);
        end
    end
    events
        tf_processor_data_ready
    end
end