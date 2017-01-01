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
        p_separations
        
        % d: calculated data
        d_time_points   double
        d_data_points   double
        d_z_values      double
        d_w_values      double
    end
    methods(Access = public)
        function this = anp_tf_processor()
            this.g_uid =            java.util.UUID.randomUUID.toString;
            this.g_wait =           false;
            this.g_busy =           false;
            this.d_z_values =       [];
            this.d_w_values =       [];
            
            tools.dbg('anp_tf_processor[constructor]:\t%s: Instance created.\n',this.g_uid);
        end
        
        function [] = init_params(this,new_params)
            this.g_wait =       true;
            if any([this.set_time_params(new_params.time_params),...
                    this.set_tf(new_params.tf_obj),...
                    this.set_angles(new_params.angles),...
                    this.set_radii(new_params.radii),...
                    this.set_weights(new_params.weights),...
                    this.set_separations(new_params.separations)])
                this.recalculate();
            end
            this.g_wait =       false;
        end
        
        function dirty = set_time_params(this,new_time_params)
            assert(~isempty(new_time_params) && isstruct(new_time_params));
            
            if isempty(this.p_time) || ~all(this.p_time.n_time_steps == new_time_params.n_time_steps)
                tools.dbg('anp_tf_processor[set_time_params]:\tNew time params.\n');
                this.p_time.n_time_steps =      new_time_params.n_time_steps;
                
                dirty =         true;
                if ~this.g_wait
                    this.recalculate();
                end
            else
                dirty =         false;
                tools.dbg('anp_tf_processor[set_time_params]:\tNothing changed.\n');
            end
        end
        
        function dirty = set_tf(this,new_tf)
            assert(~isempty(new_tf) && isa(new_tf,'tf'));
            
            if ~isequal(this.tf_obj,new_tf)
                tools.dbg('anp_tf_processor[set_tf]:\tNew transfer function\n');
                this.tf_obj =   new_tf;
                this.tf_poles = roots(new_tf.Denominator{1})';
                this.tf_zeros = roots(new_tf.Numerator{1})';
                this.tf_G =     @(z) polyval(new_tf.Numerator{1},z)./polyval(new_tf.Denominator{1},z);
                
                dirty =         true;
                if ~this.g_wait
                    this.recalculate();
                end
            else
                dirty =         false;
                tools.dbg('anp_tf_processor[set_tf]:\tNothing changed.\n');
            end
            
            this.echo_tf();
        end
        
        function dirty = set_radii(this,new_radii)
            assert(~isempty(new_radii) && isstruct(new_radii));
            
            if ~isempty(this.p_radii) && new_radii.auto_main_R
                assert(~isempty(this.p_angles) && (length(this.tf_poles) - length(this.tf_zeros) >= 0));
                possibly_new_R = anp_calc_main_R(this.tf_poles,this.tf_zeros,this.p_angles.min_angle_contribution_at_R);
            else
                possibly_new_R = NaN;
            end
            
            
            if isempty(this.p_radii) || (~isnan(possibly_new_R) && ~(this.p_radii.R == possibly_new_R)) || (isnan(possibly_new_R) && ~(this.p_radii.R == new_radii.R))
                tools.dbg('anp_tf_processor[set_radii]:\tNew radii.\n');
                this.p_radii =  new_radii;
                dirty =         true;
                if ~this.g_wait
                    this.recalculate();
                end
            else
                dirty =         false;
                tools.dbg('anp_tf_processor[set_radii]:\tNothing changed.\n');
            end
        end
        
        function dirty = set_angles(this,new_angles)
            assert(~isempty(new_angles) && isstruct(new_angles));
            
            if isempty(this.p_angles) || ~all([this.p_angles.crop_inf_transition ==         new_angles.crop_inf_transition,...
                                               this.p_angles.min_angle_contribution_at_R == new_angles.min_angle_contribution_at_R...
                                               this.p_angles.detour ==                      new_angles.detour])
                tools.dbg('anp_tf_processor[set_angles]:\tNew angles.\n');
                dirty =         true;
                this.p_angles =  new_angles;
                if ~this.g_wait
                    this.recalculate();
                end
            else
                dirty =         false;
                tools.dbg('anp_tf_processor[set_angles]:\tNothing changed.\n');
            end
        end
        
        function dirty = set_weights(this,new_weights)
            assert(~isempty(new_weights) && isstruct(new_weights));
            
            if ~isequal(this.p_weights,new_weights)
                tools.dbg('anp_tf_processor[set_weights]:\tNew weights.\n');
                dirty =             true;
                this.p_weights =    new_weights;
                if ~this.g_wait
                    this.recalculate();
                end
            else
                dirty =             false;
                tools.dbg('anp_tf_processor[set_weights]:\tNothing changed.\n');
            end
        end
        
        function dirty = set_separations(this,new_separations)
            assert(~isempty(new_separations) && isstruct(new_separations));
            
            if ~isequal(this.p_separations,new_separations)
                tools.dbg('anp_tf_processor[set_separations]:\tNew separations.\n');
                dirty =                 true;
                this.p_separations =    new_separations;
                if ~this.g_wait
                    this.recalculate();
                end
            else
                dirty =                 false;
                tools.dbg('anp_tf_processor[set_separations]:\tNothing changed.\n');
            end
        end
        
        % put the transfer function to the text output
        function [] = echo_tf(this)
            fprintf('Plotting the ');
            transfer_function = this.tf_obj %#ok<NOPRT,NASGU>
            fprintf('with\n');
            poles = this.tf_poles           %#ok<NOPRT,NASGU>
            fprintf('and\n');
            zeros = this.tf_zeros           %#ok<NOPRT,NASGU>

            if any([this.tf_obj.IODelay,this.tf_obj.InputDelay,this.tf_obj.OutputDelay])
                warning('Delay has been provided, but will be ignored, as this is currently unsupported.');
            end
        end
        
        function R = get_R(this)
            R = this.p_radii.R;
        end
        
        function data = get_time_points(this)
            data.time_points =   this.d_time_points;
            data.data_points =   this.d_data_points;
            data.time_props =    this.p_time;
        end
        
        function data = get_data(this)
            data.z_values =     this.d_z_values;
            data.w_values =     this.d_w_values;
        end
        
        function delete(this)
            tools.dbg('anp_tf_processor[delete]:\t%s: Deletion requested.\n',this.g_uid);
        end
    end
    methods(Access = private)
        function [] = recalculate(this)
            this.g_busy =       true;
            
            if this.p_radii.auto_main_R
                if isempty([this.tf_poles,this.tf_zeros])
                    this.p_radii.R = 5;
                else
                    this.p_radii.R =    anp_calc_main_R(this.tf_poles,this.tf_zeros,this.p_angles.min_angle_contribution_at_R,max(this.p_separations.pole_max,this.p_separations.zero_max));
                end
            end
            
            this.calc_time_params();
            
            args = this.prepare_z_fct_args();
            this.d_z_values = anp_in_fct_init_and_evaluate(this.d_data_points,args);
            this.d_w_values = this.tf_G(this.d_z_values);
            
            this.g_busy =   false;
        end
        
        function args = prepare_z_fct_args(this)
            args.poles =                this.tf_poles;
            args.zeros =                this.tf_zeros;
            args.radii.inf =            this.p_radii.R;
            args.angles.crop =          this.p_angles.crop_inf_transition * pi/180;
            args.angles.detour =        this.p_angles.detour * pi/180;
            args.angles.min_angle_contribution_at_R = this.p_angles.min_angle_contribution_at_R * pi/180;
            args.separation_pole_max =  this.p_separations.pole_max;
            args.separation_zero_max =  this.p_separations.zero_max;
            args.separation_margin =    this.p_separations.margin;
        end
        
        function [] = calc_time_params(this)
            tol = 100*eps;
            
            p = this.tf_poles;
            z = this.tf_zeros;
                p_neg_real =    sum(real(p) <= tol);
                p_pos_real =    sum(real(p) >= tol);
                p_pure_imag =   sum(abs(real(p)) < tol);
                z_neg_real =    sum(real(z) <= tol);
                z_pos_real =    sum(real(z) >= tol);
                z_pure_imag =   sum(abs(real(z)) < tol);
            
            total_tf_rel_deg =                  abs(p_neg_real - z_neg_real - p_pos_real + z_pos_real);
            this.p_time.oversampling_factor =   max(2,ceil(total_tf_rel_deg + this.p_weights.pole * p_pure_imag + this.p_weights.zero * z_pure_imag + floor(this.p_radii.R/10)));
            this.p_time.n_data_points =         this.p_time.n_time_steps * this.p_time.oversampling_factor;
            
            this.d_time_points =    linspace(0,1,this.p_time.n_time_steps+1);
            this.d_time_points =    this.d_time_points(1:end-1);
            this.d_data_points =    linspace(0,1,this.p_time.n_data_points);
        end
    end
end