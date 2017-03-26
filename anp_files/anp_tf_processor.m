%   class: anp_tf_processor
%   -----------------------
%   
%   This class calculates all the data that an anp_gui object then fetches
%   and subsequently plots.
%   
%   The important steps are
%   1) Receiving the main parameters and transfer function
%   2) Calculating the depending parameters like time interval size and
%      spatial oversampling
%   3) Constructing an appropriate D-Contour
%   4) Scaling the spatial resolution based on known conditions like
%      poles on the imaginary axis
%   5) Evaluating the transfer function
%   6) Providing a calling function with the precalculated data
%   
%   paper reference: p.41
classdef anp_tf_processor < handle
    
    % ---------------------------------------------------------------------
    % PRIVATE Variables
    % ---------------------------------------------------------------------
        
    properties(SetAccess = private)
        % g: general
        g_uid           % Java String,    uniqe identifier for every object instance
        s_wait          logical         % used to delay the invocation of a 'recalculate'
        s_busy          logical         % blocks any setter methods to change parameters during calculation
        
        % tf: transfer function properties
        tf_obj                          % the tf-object that we get set with 'set_tf'
        tf_poles        double          % an array of the poles that come from 'tf_obj'
        tf_zeros        double          % an array of the zeros that come from 'tf_obj'
        tf_delay        double          % if a delay was specified, we store it here
        tf_G            % function handle, is defined by 'tf_poles' and 'tf_zeros'
        
        % p: parameters
        p_time          % struct, members: 'n_time_steps','n_data_points','p_oversampling_factor',
        p_radii         % struct, members: 'auto_main_R','R'
        p_angles        % struct, members: 'crop_inf_transition','detour','min_angle_contribution'
        p_weights       % struct, members: 'pole','zero'
        p_separations   % struct, members: 'margin','pole_max','zero_max'
        
        % d: calculated data
        d_time_points   double          % 
        d_data_points   double          % 
        d_z_values      double          % data of the left plot
        d_w_values      double          % data of the right plot
    end
    
    
    
    % ---------------------------------------------------------------------
    % PUBLIC methods
    % ---------------------------------------------------------------------
    
    methods(Access = public)
        
        function this = anp_tf_processor()
            % Instantiates a tf-processor object.

            % Set a unique ID for every instance of this object
            this.g_uid =            java.util.UUID.randomUUID.toString;
            
            % Initialize the state variables
            this.s_wait =           false;
            this.s_busy =           false;
            
            tools.dbg('anp_tf_processor[constructor]:\t%s: Instance created.\n',this.g_uid);
        end
        
        function [] = init_params(this,new_params)
            % Takes all the relevant parameters and stores them in internal variables
            
            % Prevent the setter methods to invoke 'recalculate' by
            % themselves
            this.s_wait =       true;
            
            % The setter methods return a true boolean if any of their
            % arguments lead to a changed internal variable
            if any([this.set_time_params(new_params.time_params),...
                    this.set_tf(new_params.tf_obj),...
                    this.set_angles(new_params.angles),...
                    this.set_radii(new_params.radii),...
                    this.set_weights(new_params.weights),...
                    this.set_separations(new_params.separations)])
                this.recalculate();
            end
            
            % Unblock the setter methods from calling 'recalculate'
            this.s_wait =       false;
        end
        
        function dirty = set_time_params(this,new_time_params)
            % Updates the internal variables containing 'time'-information
            % Informs the caller whether something has changed
            % (dirty=true), potentially triggering a fetch of updated data
            
            assert(~isempty(new_time_params) && isstruct(new_time_params));
            
            if isempty(this.p_time) || ~all(this.p_time.n_time_steps == new_time_params.n_time_steps)
                tools.dbg('anp_tf_processor[set_time_params]:\tNew time params.\n');
                this.p_time.n_time_steps =      new_time_params.n_time_steps;
                
                dirty =         true;
                if ~this.s_wait
                    this.recalculate();
                end
            else
                dirty =         false;
                tools.dbg('anp_tf_processor[set_time_params]:\tNothing changed.\n');
            end
        end
        
        function dirty = set_tf(this,new_tf)
            % Updates the internal variables containing 'tranfer function'-information
            % Informs the caller whether something has changed
            % (dirty=true), potentially triggering a fetch of updated data
            
            assert(~isempty(new_tf) && isa(new_tf,'tf'));
            
            if ~isequal(this.tf_obj,new_tf)
                tools.dbg('anp_tf_processor[set_tf]:\tNew transfer function\n');
                this.tf_obj =   new_tf(1,1);
                this.tf_poles = roots(this.tf_obj.Denominator{1})';
                this.tf_zeros = roots(this.tf_obj.Numerator{1})';
                
                % If any sort of delay was specified, prefer IODelay over
                % Input- and OutputDelay. Calculate the effecttive delay if
                % both In- and OutputDelay are present. (According to
                % https://ch.mathworks.com/help/control/ug/time-delays-in-linear-systems.html )
                if this.tf_obj.IODelay ~= 0
                    this.tf_delay = this.tf_obj.IODelay;
                else
                    this.tf_delay = this.tf_obj.InputDelay + this.tf_obj.OutputDelay;
                end
                
                % Prepare the function that is to be evaluated later
                this.tf_G =     @(z) polyval(this.tf_obj.Numerator{1},z)./polyval(this.tf_obj.Denominator{1},z);
                
                dirty =         true;
                if ~this.s_wait
                    this.recalculate();
                end
            else
                dirty =         false;
                tools.dbg('anp_tf_processor[set_tf]:\tNothing changed.\n');
            end
            
            this.echo_tf();
        end
        
        function dirty = set_radii(this,new_radii)
            % Updates the internal variables containing information about radii of the D-contour
            
            assert(~isempty(new_radii) && isstruct(new_radii));
            
            if ~isempty(this.p_radii) && new_radii.auto_main_R
                assert(~isempty(this.p_angles) && (length(this.tf_poles) - length(this.tf_zeros) >= 0));
                possibly_new_R = anp_calc_main_R(this.tf_poles,this.tf_zeros,this.p_angles.min_angle_contribution_at_R,max(this.p_separations.pole_max,this.p_separations.zero_max));
            else
                possibly_new_R = NaN;
            end
            
            
            if isempty(this.p_radii) || (~isnan(possibly_new_R) && ~(this.p_radii.R == possibly_new_R)) || (isnan(possibly_new_R) && ~(this.p_radii.R == new_radii.R))
                tools.dbg('anp_tf_processor[set_radii]:\tNew radii.\n');
                this.p_radii =  new_radii;
                dirty =         true;
                if ~this.s_wait
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
                if ~this.s_wait
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
                if ~this.s_wait
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
                if ~this.s_wait
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
            transfer_function = this.tf_obj(1,1)        %#ok<NOPRT,NASGU>
            fprintf('with\n');
            poles = this.tf_poles                       %#ok<NOPRT,NASGU>
            fprintf('and\n');
            zeros = this.tf_zeros                       %#ok<NOPRT,NASGU>
            
            if this.tf_delay ~= 0
                fprintf('and\n');
                delay = this.tf_delay                   %#ok<NOPRT,NASGU>
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
            this.s_busy =       true;
            
            if this.p_radii.auto_main_R
                if isempty([this.tf_poles,this.tf_zeros])
                    this.p_radii.R = 5;
                else
                    this.p_radii.R =    anp_calc_main_R(this.tf_poles,this.tf_zeros,this.p_angles.min_angle_contribution_at_R,max(this.p_separations.pole_max,this.p_separations.zero_max),this.tf_delay);
                end
            end
            
            this.calc_time_params();
            
            args = this.prepare_z_fct_args();
            this.d_z_values = anp_in_fct_init_and_evaluate(this.d_data_points,args);
            
            if this.tf_delay ~= 0
                this.d_w_values = exp(imag(this.d_z_values) * this.tf_delay * 1i) .* this.tf_G(this.d_z_values);
            else
                this.d_w_values = this.tf_G(this.d_z_values);
            end
            
            this.s_busy =   false;
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
            
            radius_correction =     floor(this.p_radii.R/10);
            
            pz = [p,z];
            pz_not_on_im_axis = pz(abs(real(pz)) >= 100*eps);
            %pz_min_abs = min(abs([1,pz_not_on_im_axis]));
            
            n_small_poles = sum(abs(p) < 1 & abs(real(p)) >= 100*eps);
            n_small_zeros = sum(abs(z) < 1 & abs(real(z)) >= 100*eps);
            dist = mean(abs(pz_not_on_im_axis(abs(pz_not_on_im_axis) < 0.1)));
            small_pz_correction =   (1+radius_correction) * max(0, n_small_poles - n_small_zeros) * (-5*log(dist));
            %old_pz_correction =     -floor(20*log10(pz_min_abs));
            
            im_pole_correction =    this.p_weights.pole * p_pure_imag;
            im_zero_correction =    this.p_weights.zero * z_pure_imag;
            
            this.p_time.oversampling_factor =   max(5,ceil(total_tf_rel_deg + im_pole_correction + im_zero_correction + radius_correction + small_pz_correction + this.tf_delay));
            
            % trade off oversampling for time steps (animation speed will
            % be lower) if oversampling is high
            tradeoff = fix(1 + this.p_time.oversampling_factor/100);
            this.p_time.oversampling_factor = fix(this.p_time.oversampling_factor / tradeoff);
            this.p_time.n_time_steps = this.p_time.n_time_steps * tradeoff;
            this.p_time.n_data_points =         this.p_time.n_time_steps * this.p_time.oversampling_factor;
            
            this.d_time_points =    linspace(0,1,this.p_time.n_time_steps+1);
            this.d_time_points =    this.d_time_points(1:end-1);
            this.d_data_points =    linspace(0,1,this.p_time.n_data_points);
        end
    end
end