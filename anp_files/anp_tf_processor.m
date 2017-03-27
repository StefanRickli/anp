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
%   Parametrization terminology:
%   We calculate for every 'data-point' (in [0,1]) the corresponding values
%   in the z- and w plane.
%   As we require to maintain a certain maximal spacing between w-points,
%   it's no good idea to iterate through all of them in the animation.
%   Instead we only take every 'oversampling_factor'-th point. This is the
%   'time-point' parametrization.
%   
%   We use the following mapping:
%
%                                   nonlinear mon. increasing map                                  D-contour shapes            G(z)
%   data_points (t) = linspace(0,1) -----------------------------> q [0,total_shape_circumference] ----------------> z complex ----> w complex
%   
%   We need the pre-mapping from t to q in order to compensate the
%   highly nonlinear relationship between the animation parametetrization
%   and the speed the arrow moves with in the Nyquist plot (w-plane)
%   For example: while we take a detour around a pole on the imaginary
%                axis, the distance that we take on the D-contour is tiny,
%                but in the Nyquist plot the curve has a lot of movement.
%                If we were to simply assign the same evenly spaced points
%                for the detour as for the preceeding straight line on the
%                imaginary axis, we would get very poor spatial resolution
%                in the Nyquist plot.
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
            % 
            % As all the called functions have side effects that we rely on
            % we need that all of them really run. Hence we can't use a
            % syntax like A || B || C etc because of lazy computing.
            if any([this.set_time_params(new_params.time_params),...
                    this.set_tf(new_params.tf_obj,new_params.tf_delay),...
                    this.set_angles(new_params.angles),...
                    this.set_separations(new_params.separations),...
                    this.set_radii(new_params.radii),...
                    this.set_weights(new_params.weights),])
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
        
        function dirty = set_tf(this,new_tf,delay)
            % Updates the internal variables containing 'tranfer function'-information
            % Informs the caller whether something has changed
            % (dirty=true), potentially triggering a 'recalculate' and/or
            % a fetch of updated data
            
            assert(~isempty(new_tf) && isa(new_tf,'tf'));
            
            if ~isequal(this.tf_obj,new_tf)
                tools.dbg('anp_tf_processor[set_tf]:\tNew transfer function\n');
                this.tf_obj =   new_tf;
                this.tf_poles = roots(this.tf_obj.Denominator{1})';
                this.tf_zeros = roots(this.tf_obj.Numerator{1})';
                
                this.tf_delay = delay;
                
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
            % Informs the caller whether something has changed
            % (dirty=true), potentially triggering a 'recalculate' and/or
            % a fetch of updated data
            
            assert(~isempty(new_radii) && isstruct(new_radii));
            
            % auto_main_R == true whenever a user hasn't set a manual value
            % for R
            if new_radii.auto_main_R
                % 'anp_calc_main_R' needs some angles and the separations
                % to be set in order to work properly. Also we want a
                % proper transfer function.
                assert(~isempty(this.p_angles) && ~isempty(this.p_separations) && (length(this.tf_poles) - length(this.tf_zeros) >= 0));
                
                % Check whether the transfer function is a constant
                if isempty([this.tf_poles,this.tf_zeros])
                    this.p_radii.R = 5;
                else
                    new_radii.R =    anp_calc_main_R(this.tf_poles,this.tf_zeros,this.p_angles.min_angle_contribution_at_R,max(this.p_separations.pole_max,this.p_separations.zero_max),this.tf_delay);
                end
            end
            
            if isempty(this.p_radii) || this.p_radii.R ~= new_radii.R
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
            % Updates the internal variables containing information about angles of the D-contour
            % Informs the caller whether something has changed
            % (dirty=true), potentially triggering a 'recalculate' and/or
            % a fetch of updated data
            
            assert(~isempty(new_angles) && isstruct(new_angles));
            
            if isempty(this.p_angles) || ~all([this.p_angles.crop_inf_transition ==         new_angles.crop_inf_transition,...
                                               this.p_angles.min_angle_contribution_at_R == new_angles.min_angle_contribution_at_R...
                                               this.p_angles.detour ==                      new_angles.detour])
                tools.dbg('anp_tf_processor[set_angles]:\tNew angles.\n');
                this.p_angles = new_angles;
                dirty =         true;
                if ~this.s_wait
                    this.recalculate();
                end
            else
                dirty =         false;
                tools.dbg('anp_tf_processor[set_angles]:\tNothing changed.\n');
            end
        end
        
        function dirty = set_weights(this,new_weights)
            % Updates the internal variables containing information about weigths that influence the spatial resolution of the D-contour
            % Informs the caller whether something has changed
            % (dirty=true), potentially triggering a 'recalculate' and/or
            % a fetch of updated data
            
            assert(~isempty(new_weights) && isstruct(new_weights));
            
            if ~isequal(this.p_weights,new_weights)
                tools.dbg('anp_tf_processor[set_weights]:\tNew weights.\n');
                this.p_weights =    new_weights;
                dirty =             true;
                if ~this.s_wait
                    this.recalculate();
                end
            else
                dirty =             false;
                tools.dbg('anp_tf_processor[set_weights]:\tNothing changed.\n');
            end
        end
        
        function dirty = set_separations(this,new_separations)
            % Updates the internal variables containing information about separations of detours on the imaginary axis of the D-contour
            % Informs the caller whether something has changed
            % (dirty=true), potentially triggering a 'recalculate' and/or
            % a fetch of updated data
            
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
        
        function [] = echo_tf(this)
            % Puts information about the transfer function to the text output
            
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
            % Access method to get R
            
            R = this.p_radii.R;
        end
        
        function data = get_time_points(this)
            % Access method to get the parametrization
            
            data.time_points =   this.d_time_points;
            data.data_points =   this.d_data_points;
            data.time_props =    this.p_time;
        end
        
        function data = get_data(this)
            % Access method to get the data for the (left) z-plot and (right) w-plot
            
            data.z_values =     this.d_z_values;
            data.w_values =     this.d_w_values;
        end
        
        function delete(this)
            % Destructor method
            
            tools.dbg('anp_tf_processor[delete]:\t%s: Deletion requested.\n',this.g_uid);
        end
    end
    
    % ---------------------------------------------------------------------
    % PRIVATE methods
    % ---------------------------------------------------------------------
    
    methods(Access = private)
        function [] = recalculate(this)
            % Calculates the w- and z-data
            % 
            % This is the heart of the tf_processor object!
            % PRE: all the necessary parameters by the user have been set.
            % 
            % TODO: potential bug: as the transfer function can be set
            % independently from R, R possibly doesn't get updated before
            % recalculation
            
            this.s_busy =       true;
            
            % Based on the parameters (transfer function, R, etc),
            % calculate how many evaluations of G(s) we need and store the
            % data-point parametrization in the internal variable
            this.calc_parametrization();
            
            % Put all the arguments together that we then pass to the
            % function that maps the data-points from [0,1] to the
            % D-contour in the z-plot
            args = this.prepare_z_fct_args();
            this.d_z_values = anp_d_contour_init_and_evaluate(this.d_data_points,args);
            
            if this.tf_delay ~= 0
                this.d_w_values = exp(-imag(this.d_z_values) * this.tf_delay * 1i) .* this.tf_G(this.d_z_values);
            else
                this.d_w_values = this.tf_G(this.d_z_values);
            end
            
            this.s_busy =   false;
        end
        
        function args = prepare_z_fct_args(this)
            % Puts together all the arguments for the function that maps [0,1] to the D-contour
            
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
        
        function [] = calc_parametrization(this)
            % Decides on the distribution of time-steps
            % 
            % The GUI will plot D-contour(data-points) and iterate through
            % the time-points for the arrow.
            
            tol = 100*eps;
            
            % Count the # of poles/zeros with zero/positive/negative real
            % part
            p = this.tf_poles;
            z = this.tf_zeros;
                p_neg_real =    sum(real(p) <= tol);
                p_pos_real =    sum(real(p) >= tol);
                p_pure_imag =   sum(abs(real(p)) < tol);
                z_neg_real =    sum(real(z) <= tol);
                z_pos_real =    sum(real(z) >= tol);
                z_pure_imag =   sum(abs(real(z)) < tol);
            
            % Calc the effective relative degree of the transfer function
            % polynomial. This has an effect on how fast the phase changes
            % in the orgin during the half-circle of the D-contour.
            total_tf_rel_deg =      abs(p_neg_real - z_neg_real - p_pos_real + z_pos_real);
            
            % If we go to high frequencies (that is a large radius of the
            % half-circle, R), we loose spatial resolution at low
            % frequencies. Compensate this.
            radius_correction =     floor(this.p_radii.R/5);
            
            % Poles that are near the imaginary axis will lead to huge
            % movement in the Nyquist plot when the w-arrow gets near them.
            % Zeros that are near those poles will compensate their effect.
            % So we only compensate their relative degree.
            % 
            % TODO: We could do better by analyzing the exact distribution
            % of those p/z with small real part, but we stick with this
            % implementation at the moment.
            pz = [p,z];
            pz_not_on_im_axis =     pz(abs(real(pz)) >= 100*eps);
            n_small_poles =         sum(abs(p) < 1 & abs(real(p)) >= 100*eps);
            n_small_zeros =         sum(abs(z) < 1 & abs(real(z)) >= 100*eps);
            magnitude =             mean(abs(pz_not_on_im_axis(abs(pz_not_on_im_axis) < 0.1)));
            small_pz_correction =   (1+radius_correction) * max(0, n_small_poles - n_small_zeros) * (-5*log(magnitude));
            
            % Poles and zeros that are exactly on the imaginary axis lead
            % to the creation of a detour in the D-contour.
            % Purely imaginary poles cause large movement in the Nyquist
            % plot during the detour.
            % Purely imaginary zeros cause large movement right before and
            % after the detour as the mapped points are forced to go to the
            % origin from wherever they were.
            im_pole_correction =    this.p_weights.pole * p_pure_imag;
            im_zero_correction =    this.p_weights.zero * z_pure_imag;
            
            % Put together all the corrections and calculate a factor that
            % relates the # of animation steps (arrow movement) with the
            % spatial resolution of the plot.
            this.p_time.oversampling_factor =   max(5,ceil(total_tf_rel_deg + im_pole_correction + im_zero_correction + radius_correction + small_pz_correction + this.tf_delay));
            
            % Trade off oversampling for time steps (animation speed will
            % be lower) if oversampling is high
            % 
            % TODO: check the rounding procedure, why already fix tradeoff?
            tradeoff =                          fix(1 + this.p_time.oversampling_factor/100);
            this.p_time.oversampling_factor =   fix(this.p_time.oversampling_factor / tradeoff);
            this.p_time.n_time_steps =          this.p_time.n_time_steps * tradeoff;
            this.p_time.n_data_points =         this.p_time.n_time_steps * this.p_time.oversampling_factor;
            
            % Finally generate the vectors of the parametrization
            this.d_time_points =    linspace(0,1,this.p_time.n_time_steps+1);
            this.d_time_points =    this.d_time_points(1:end-1);
            this.d_data_points =    linspace(0,1,this.p_time.n_data_points);
        end
    end
end