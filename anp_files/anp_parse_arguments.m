% TODO duration, FPS are deprecated
% 2017-02-24: reuse them later

function checked_args = anp_parse_arguments(varargin)
    
    %% Globals
    % Define the allowed single keywords here in 'anp_keywords'
    global anp_arg_types anp_keywords;
    anp_arg_types = '';
    anp_keywords = {'return_handle',false;'trigger_step',false;'cleanup_after_error',false};
    
    %% Argument syntax checking
    ip = inputParser;
    
    default_arg1            = NaN;
    default_arg2            = NaN;
    default_arg3            = NaN;
    default_arg4            = NaN;
    default_arg5            = NaN;
    
    default_R               = NaN;
    default_duration        = 20;       % length of animation in movie player, affects spatial resolution, [s]
    default_FPS             = 6;        % speed of animation in movie player, affects spatial resolution, [frames/s]
    default_trail_length    = 0.15;     % how long should the trail be?, [0<trail_length<1]
    default_z_plot_x0       = [0;0];    % center of left plot
    default_z_plot_dims     = [5;5];    % width and heigth of left plot
    default_w_plot_x0       = [0;0];    % center of right plot
    default_w_plot_dims     = [2;2];    % width and heigth of right plot
    default_plot_size       = 500;      % [pixel]
    default_filename        = 'animated_nyquist_plot_export';
    
    addOptional(ip,'arg1',          default_arg1,           @checkOptArg);
    addOptional(ip,'arg2',          default_arg2,           @checkOptArg);
    addOptional(ip,'arg3',          default_arg3,           @checkOptArg);
    addOptional(ip,'arg4',          default_arg4,           @checkOptArg);
    addOptional(ip,'arg5',          default_arg5,           @checkOptArg);
    
    addParameter(ip,'R',            default_R,              @isScalarNumberPositive);
    addParameter(ip,'duration',     default_duration,       @isScalarNumber);
    addParameter(ip,'FPS',          default_FPS,            @isScalarInteger);
    addParameter(ip,'trail_length', default_trail_length,   @checkTrail);
    addParameter(ip,'left_x0',      default_z_plot_x0,      @isVectorNumber);
    addParameter(ip,'left_dims',    default_z_plot_dims,    @isVectorNumber);
    addParameter(ip,'right_x0',     default_w_plot_x0,      @isVectorNumber);
    addParameter(ip,'right_dims',   default_w_plot_dims,    @isVectorNumber);
    addParameter(ip,'plot_size',    default_plot_size,      @isScalarInteger);
    addParameter(ip,'filename',     default_filename,       @ischar);
    
    % do the first round of syntax checking, the parsers also build a
    % string, 'anp_arg_types', according to the types of arguments encountered
    parse(ip,varargin{:});
    
    % for the optional arguments we require
    % - either (neither a transfer function nor poles/zeros), followed by zero, one or two anp_keywords
    % - OR a transfer function, followed by zero, one or two anp_keywords
    % - OR two vectors containing the poles and zeros, followed by zero, one or two anp_keywords
    if isempty(regexp(anp_arg_types,'^(t|v{2})?k{0,3}$','emptymatch'))
        error('Check type and order of optional arguments.');
    end
    
    %% Argument extraction

    % ---------------------------------------------------------------------
    % Transfer function parameters
    % ---------------------------------------------------------------------
    
    % unused yet
    % parse the poles, zeros and construct the tf-object if necessary
    if ~isempty(regexp(anp_arg_types,'^k{0,3}$','emptymatch','once'))
        % if neither a [t]ransfer function nor [v]ectors of poles and
        % zeros have been specified: use some standard values for
        % demonstration purposes
        
        fprintf('\nNo transfer function specified. Showing default example.\n\n');
        fprintf('Usage:\n');
        fprintf('animated_nyquist_plot(tf)\t\t\t\twhere ''tf'' is a transfer function object by Matlab command ''tf''\n');
        fprintf('animated_nyquist_plot([zeros],[poles])\twhere ''zeros'' and ''poles'' are 2D row vectors with roots of [num] and [denum]\n');
        fprintf('animated_nyquist_plot(__,Name,Value)\tsee documentation\n\n');
        
        tf_poles =                  [-3,-2,-1+1i,-1-1i];
        tf_zeros =                  -0.7;
        checked_args.tf_obj =       tf(poly(tf_zeros),poly(tf_poles));
        
    elseif ~isempty(regexp(anp_arg_types,'^tk{0,3}$', 'once'))
        % [t]ransfer function
        % we only take the first transfer function if an object with 
        % multiple transfer function has been provided
        checked_args.tf_obj =       ip.Results.arg1(1,1);
        
        if any(size(checked_args.tf_obj.Numerator)-[1,1])
            warning('Specified system is not SISO. Will only use first transfer function from input 1 to output 1.');
        end
        
    elseif ~isempty(regexp(anp_arg_types,'^vvk{0,3}$', 'once'))
        % two [v]ectors, containing lists of pole and zero locations
        
        % transpose to row vector if necessary
        if ~isempty(ip.Results.arg2) && length(ip.Results.arg2(:,1)) > 1
            tf_poles = ip.Results.arg2.';
        else
            tf_poles = ip.Results.arg2;
        end
        
        if ~isempty(ip.Results.arg1) && length(ip.Results.arg1(:,1)) > 1
            tf_zeros = ip.Results.arg1.';
        else
            tf_zeros = ip.Results.arg1;
        end
        
        checked_args.tf_obj =       tf(poly(tf_zeros),poly(tf_poles));
    else
        error('Error: oops, we shouldn''t be here... sorry about that. Please send me an email about this and provide me with the input arguments you used.');
    end
    
    % If any sort of delay was specified, prefer IODelay over
    % Input- and OutputDelay. Calculate the effective delay if
    % both In- and OutputDelay are present. (According to
    % https://ch.mathworks.com/help/control/ug/time-delays-in-linear-systems.html )
    if checked_args.tf_obj.IODelay ~= 0
        checked_args.tf_delay = checked_args.tf_obj.IODelay;
    else
        checked_args.tf_delay = checked_args.tf_obj.InputDelay + checked_args.tf_obj.OutputDelay;
    end
    
    % ---------------------------------------------------------------------
    % Plot parameters
    % ---------------------------------------------------------------------
    checked_args.trail_length =         ip.Results.trail_length;
    
    % parse the parameters for the z-plot (left)
    % the GUI will choose the parameters automatically
    % if z_plot_auto_lims = true
    if ~isequal(ip.Results.left_x0,default_z_plot_x0) || ~isequal(ip.Results.left_dims,default_z_plot_dims)
        checked_args.z_plot_auto_lims = false;
        
        if ~isempty(ip.Results.left_x0) && length(ip.Results.left_x0(1,:)) > 1
            % argument has been given but as column vector, use transposed
            checked_args.z_plot_x0 =    ip.Results.left_x0';
        else
            % either 'ip.Results.left_x0' is initialized with the default
            % value 'default_z_plot_x0' if no argument has been given
            % OR it contains the argument by the user
            checked_args.z_plot_x0 =    ip.Results.left_x0;
        end
        
        if ~isempty(ip.Results.left_dims) && length(ip.Results.left_dims(1,:)) > 1
            checked_args.z_plot_dims =  ip.Results.left_dims';
        else
            checked_args.z_plot_dims =  ip.Results.left_dims;
        end
    else
        checked_args.z_plot_auto_lims = true;
        checked_args.z_plot_x0 =        [];
        checked_args.z_plot_dims =      [];
    end
    
    % parse the parameters for the w-plot (right), see documentation for
    % z-plot
    if ~isequal(ip.Results.right_x0,default_w_plot_x0) || ~isequal(ip.Results.right_dims,default_w_plot_dims)
        checked_args.w_plot_auto_lims = false;
        
        if ~isempty(ip.Results.right_x0) && length(ip.Results.right_x0(1,:)) > 1
            checked_args.w_plot_x0 =   ip.Results.right_x0';
        else
            checked_args.w_plot_x0 =   ip.Results.right_x0;
        end
        
        if ~isempty(ip.Results.right_dims) && length(ip.Results.right_dims(1,:)) > 1        
            checked_args.w_plot_dims = ip.Results.right_dims';
        else
            checked_args.w_plot_dims = ip.Results.right_dims;
        end
    else
        checked_args.w_plot_auto_lims = true;
        checked_args.w_plot_x0 =        [];
        checked_args.w_plot_dims =      [];
    end
    
    % ---------------------------------------------------------------------
    % Figure parameters
    % ---------------------------------------------------------------------
    checked_args.plot_size =            ip.Results.plot_size;
    checked_args.border =               20;                 % pixel   
    
    % ---------------------------------------------------------------------
    % Input-function parameters
    % ---------------------------------------------------------------------
    if isnan(ip.Results.R)
        checked_args.radii.auto_main_R = true;
    else
        checked_args.radii.auto_main_R = false;
    end
    checked_args.radii.R =              ip.Results.R;
    checked_args.angles.crop_inf_transition =           7;  % [°], roundoff-parameter of the D-curve, usually not necessary to change
    checked_args.angles.min_angle_contribution_at_R =   65; % [°]
    checked_args.angles.detour =        45;                 % [°]
    checked_args.halfsecants.pole_max = 1/4;                % [1] absolute
    checked_args.halfsecants.zero_max = 1/8;                % [1] absolute
    checked_args.halfsecants.margin =   0.3;                % how much free space between the nearest neighboring poles/zeros? ==> avoid that the nearest pole/zero-detours could have no straight part between them. p.31
    checked_args.weight.per_pole =      2;                  % if a pole is on the Im-axis, how many spatial points should it be given in the nyquist plot?
    checked_args.weight.per_zero =      1/3;                % same as above
    checked_args.weight.of_axis_parts = 9;                  % 
    checked_args.weight.of_crop_inf =   3;                  % 
    
    % ---------------------------------------------------------------------
    % Time and spatial resolution parameters
    % ---------------------------------------------------------------------
    checked_args.time_params.n_time_steps =         120;    % a) how many steps in the animation? b) base value for # of data points
    checked_args.time_params.resolution_factor =    3;      % TODO deprecated, unused: remove in future version
    
    % ---------------------------------------------------------------------
    % Keywords
    % ---------------------------------------------------------------------
    checked_args.return_handle =    anp_keywords{1,2};      % Causes anp_main to return the handle to the GUI and tf_processor object
    checked_args.trigger_step =     anp_keywords{2,2};      % Causes anp_main to trigger one call to the step callback function in anp_gui
    checked_args.cleanup_after_error = anp_keywords{3,2};   % Causes anp_main to close all windows and delete its objects.
    
    % ---------------------------------------------------------------------
    % Video export parameters, DEPRECATED
    % ---------------------------------------------------------------------
    checked_args.filename =     ip.Results.filename;        % TODO deprecated, unused
    
end

%% Syntax checking functions for Matlabs input parser
function res = checkOptArg(x)
    global anp_arg_types
    
    % simply try to find out the type of argument that has been provided
    % do a check for sanity later
    n = length(anp_arg_types);
    if isa(x,'tf')
        anp_arg_types(n+1) = 't';                 % t = transfer function
    elseif isempty(x) || (isvector(x) && isnumeric(x))
        anp_arg_types(n+1) = 'v';                 % v = vector
    elseif isKeyword(x)
        anp_arg_types(n+1) = 'k';                 % k = keyword
    else
        res = false;
        return;
    end
    
    res = true;
    return;
end

function res = isKeyword(x)
    % Checks whether the provided keyword is a valid one. If it is, it remembers that the keyword was set by setting the corresponding boolean to True.
    
    global anp_keywords;
    
    try
        res = ismember(lower(x),lower({anp_keywords{:,1}}));
    catch
        error('Wrong argument type.');
    end
    
    for ii = 1:length({anp_keywords{:,1}})
        anp_keywords{ii,2} = anp_keywords{ii,2} || isequal(x,anp_keywords{ii,1});
    end
end

function res = isScalarInteger(x)
    res = isscalar(x) && ~mod(x,1);
end

function res = isScalarNumber(x)
    res = isscalar(x) && isnumeric(x);
end

function res = checkTrail(x)
    res = isscalar(x) && isnumeric(x) && x > 0 && x < 1;
end

function res = isScalarNumberPositive(x)
    res = isscalar(x) && isnumeric(x) && x > 0;
end

function res = isVectorNumber(x)
    res = (isequal(size(x),[1,2]) || isequal(size(x),[2,1])) && isnumeric(x);
end
