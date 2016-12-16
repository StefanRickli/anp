% -------------------------------------------------------------------------
%   Animated Nyquist Plot
%   *********************
%
%   Author: Stefan Rickli, ricklis [at] student.ethz.ch
%           https://blogs.ethz.ch/ricklis
%
%   Version 3.2, 2016-11-23
%   
%   
%   What does it do?:
%   -----------------
%   This script calculates for a given set of zeros and poles a suitable
%   D-curve, plots it together with the zeros and poles and draws the
%   corresponding nyquist curve in another plot.
%   For every time step one is shown the contribution of each zero and pole
%   to the magnitude and phase of the nyquist curve's values.
%   After calculation, the animation is loaded into an interactive player
%   where one can further investigate each step.
%   
%   Syntax:
%   -------
%   animated_nyquist_plot(tf)
%   animated_nyquist_plot([zeros],[poles])
%   animated_nyquist_plot(__,Name,Value)
%   
%   tf: An object of type 'transfer function', created by the command 'tf'.
%       Note that only the first transfer function is considered, should
%       there be multiple.
%   
%   [zeros],[poles]:    Row arrays with the roots of the Numerator and
%                       Denominator of the transfer function
%   
%   Name-Value pairs:
%   Per default the function manages its parameters itself. But you can
%   specify certain aspects manually:
%   
%   'R':    Sets the radius (that is also about the maximum omega) of the 
%           D-curve manually. Per default it choses its size twice the
%           absolute value of the largest pole or zero. See limitation (1).
%   
%   'plot_size':    Specifies the dimensions of each of the plots in
%                   pixels.
%                   Default size is 500px.
%   
%   'duration': How long should the animation in the video player take.
%               Lower values:   shorter computation time
%                               less spatial resolution
%               Unit: seconds
%   
%   'FPS':      Determines the FPS (frames per second) of the animation.
%               Lower values:   shorter computation time
%                               less spatial resolution
%   
%   'trail-length': Determines the fraction of the plot that is colored
%                   orange. Acceptable values are 0 <= length < 1.
%   
%   'left_x0':  Manual value for the center of the left subplot.
%               Provide a 2D row vector, e.g. [x0,y0].
%               If no width and height are specified by 'left_dims', we
%               assume [5,5] per default.
%   
%   'left_dims':    Manual value for the width and height of the left
%                   subplot.
%                   Provide a 2D row vector, e.g. [width,height].
%                   If no center is specified by 'left_x0', we assume
%                   [0,0] per default.
%   
%   'right_x0': Same functionality as 'left_x0', but for the right subplot.
%               Default 'right_dims' if not provided: [2,2]
%   
%   'right_dims':   Same functionality as 'left_dims', but for the right
%                   subplot.
%                   Default 'right_x0' if not provided: [0,0]
%   
%   
%   Known Limitations:
%   ------------------
%   1)  Although the radius of the D-curve adapts itself to the magnitude
%       of the largest pole or zero, it doesn't correct large poles or
%       zeros that are near the real axis (yet). This causes the maximum
%       angle contribution to be less and less the bigger the real part of
%       said pole or zero is. 
%   2)  The script can't handle poles and zeros that lie exactly on the
%       imaginary axis. Please place those points slightly to the left
%       (or right).
%   3)  At the moment the function can't handle scaling factors in front of
%       the transfer function, e.g. the k in G(s)=k*(s+1)/(s+2)
%
%
%   Known Issues:
%   -------------
%
%   - axis numbers sometimes disappear with certain combinations of
%     zeros/poles
%
% -------------------------------------------------------------------------

function [] = animated_nyquist_plot_new(varargin)
    ip = inputParser;
    
    global arg_types keywords;
    arg_types = '';
    keywords = {'noplayer',false;'export',false};
    
    default_arg1            = NaN;
    default_arg2            = NaN;
    default_arg3            = NaN;
    default_arg4            = NaN;
    
    default_R               = NaN;
    default_duration        = 20;               % length of animation in movie player, affects spatial resolution, [s]
    default_FPS             = 18;               % speed of animation in movie player, affects spatial resolution, [frames/s]
    default_trail_length    = 0.15;             % how long should the trail be?, [0<trail_length<1]
    default_left_x0         = [0;0];            % center of left plot
    default_left_dims       = [5;5];            % width and heigth of left plot
    default_right_x0        = [0;0];            % center of right plot
    default_right_dims      = [2;2];            % width and heigth of right plot
    default_plot_size       = 500;              % pixel
    default_filename        = 'animated_nyquist_plot_export';
    
    addOptional(ip,'arg1',          default_arg1,           @checkOptArg);
    addOptional(ip,'arg2',          default_arg2,           @checkOptArg);
    addOptional(ip,'arg3',          default_arg3,           @checkOptArg);
    addOptional(ip,'arg4',          default_arg4,           @checkOptArg);
    
    addParameter(ip,'R',            default_R,              @isScalarNumberPositive);
    addParameter(ip,'duration',     default_duration,       @isScalarNumber);
    addParameter(ip,'FPS',          default_FPS,            @isScalarInteger);
    addParameter(ip,'trail_length', default_trail_length,   @checkTrail);
    addParameter(ip,'left_x0',      default_left_x0,        @isVectorNumber);
    addParameter(ip,'left_dims',    default_left_dims,      @isVectorNumber);
    addParameter(ip,'right_x0',     default_right_x0,       @isVectorNumber);
    addParameter(ip,'right_dims',   default_right_dims,     @isVectorNumber);
    addParameter(ip,'plot_size',    default_plot_size,      @isScalarInteger);
    addParameter(ip,'filename',     default_filename,       @ischar);
    
    parse(ip,varargin{:});
    
    % for the optional arguments we require
    % either (neither a transfer function nor poles/zeros), followed by zero, one or two keywords
    % OR a transfer function, followed by zero, one or two keywords
    % OR two vectors containing the poles and zeros, followed by zero, one or two keywords
    if isempty(regexp(arg_types,'^(t|v{2})?k{0,2}$','emptymatch'))
        error('Check type and order of optional arguments.');
    end
    
    if ~isempty(regexp(arg_types,'^k{0,2}$','emptymatch','once'))
        % if no poles and zeros have been specified
        fprintf('\nNo transfer function specified. Showing default example.\n\n');
        fprintf('Usage:\n');
        fprintf('animated_nyquist_plot(tf)\t\t\t\twhere ''tf'' is a transfer function object by Matlab command ''tf''\n');
        fprintf('animated_nyquist_plot([zeros],[poles])\twhere ''zeros'' and ''poles'' are 2D row vectors with roots of [num] and [denum]\n');
        fprintf('animated_nyquist_plot(__,Name,Value)\tsee documentation\n');
        
        functionParams.tf_zeros = [-0.7];
        functionParams.tf_poles = [-3,-2,-1+1i,-1-1i];
        
    elseif ~isempty(regexp(arg_types,'^tk{0,2}$', 'once'))
        % we only take the first transfer function if there multiple have
        % been provided
        functionParams.tf_zeros = roots(ip.Results.arg1.Numerator{1})';         % zeros of transfer function
        functionParams.tf_poles = roots(ip.Results.arg1.Denominator{1})';       % poles of transfer function
        
    elseif ~isempty(regexp(arg_types,'^vvk{0,2}$', 'once'))
        if ~isempty(ip.Results.arg1) && length(ip.Results.arg1(:,1)) > 1
            functionParams.tf_zeros = ip.Results.arg1';
        else
            functionParams.tf_zeros = ip.Results.arg1;
        end
        
        if ~isempty(ip.Results.arg2) && length(ip.Results.arg2(:,1)) > 1
            functionParams.tf_poles = ip.Results.arg2';
        else
            functionParams.tf_poles = ip.Results.arg2
        end
    else
        error('Error: oops, we shouldn''t be here... sorry about that. Please send me an email about this and provide me with the input arguments you used.');
    end
    
    if ~isequal(ip.Results.left_x0,[0;0]) || ~isequal(ip.Results.left_dims,[5;5])
        functionParams.in_auto_size = false;        % values below are only active if auto = false

        if ~isempty(ip.Results.left_x0) && length(ip.Results.left_x0(1,:)) > 1
            functionParams.in_x0 = ip.Results.left_x0';
        else
            functionParams.in_x0 = ip.Results.left_x0;
        end
        
        if ~isempty(ip.Results.left_dims) && length(ip.Results.left_dims(1,:)) > 1
            functionParams.in_dims = ip.Results.left_dims';
        else
            functionParams.in_dims = ip.Results.left_dims;
        end
    else
        functionParams.in_auto_size = true;
    end
    
    if ~isequal(ip.Results.right_x0,[0;0]) || ~isequal(ip.Results.right_dims,[2;2])
        functionParams.out_auto_size = false;       % values below are only active if auto = false
        
        if ~isempty(ip.Results.right_x0) && length(ip.Results.right_x0(1,:)) > 1
            functionParams.out_x0 = ip.Results.right_x0';
        else
            functionParams.out_x0 = ip.Results.right_x0;
        end
        
        if ~isempty(ip.Results.right_dims) && length(ip.Results.right_dims(1,:)) > 1        
            functionParams.out_dims = ip.Results.right_dims';
        else
            functionParams.out_dims = ip.Results.right_dims;
        end
    else
        functionParams.out_auto_size = true;
    end
    
    functionParams.duration = ip.Results.duration;
    functionParams.FPS = ip.Results.FPS;
    functionParams.trail_length = ip.Results.trail_length;
    functionParams.plotSize = ip.Results.plot_size;
    functionParams.border = 20;                    % pixel   
    functionParams.noPlayer = keywords{1,2};
    functionParams.export = keywords{2,2};
    functionParams.filename = ip.Results.filename;
    functionParams.R = ip.Results.R;
    functionParams.min_angle_contribution_at_R = 65;
    
    
    main(functionParams);
end

function res = checkOptArg(x)
    global arg_types
    
    % simply tries to find out the type of argument that has been provided
    % do a check for sanity later
    n = length(arg_types);
    if isa(x,'tf')
        arg_types(n+1) = 't';                 % t = transfer function
    elseif isempty(x) || (isvector(x) && isnumeric(x))
        arg_types(n+1) = 'v';                 % v = vector
    elseif isKeyword(x)
        arg_types(n+1) = 'k';                 % k = keyword
    else
        res = false;
        return;
    end
    
    res = true;
    return;
end

function res = isKeyword(x)
    global keywords;
    
    try
        res = ismember(lower(x),lower({keywords{:,1}}));
    catch
        error('Wrong argument type.');
    end
    
    for ii = 1:length({keywords{:,1}})
        keywords{ii,2} = keywords{ii,2} || isequal(x,keywords{ii,1});
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

%% main function

function [] = main(animationParams)
    %% debug
    animationParams.t_min = 0;
    animationParams.t_max = 1;
        
    % put the transfer function to the text output
    transfer_function = tf(poly(animationParams.tf_zeros),poly(animationParams.tf_poles))
    fprintf('with\n');
    poles = animationParams.tf_poles
    fprintf('and\n');
    zeros = animationParams.tf_zeros
    
    %keyboard;
    
    %% initialization
    animationParams.phi = 7;                                    % [°], roundoff-parameter of the D-curve, usually not necessary to change
    if isnan(animationParams.R)
        pz_non_only_imaginary = [poles,zeros];
        R1 = imag(pz_non_only_imaginary) - real(pz_non_only_imaginary) * tan(deg2rad(animationParams.min_angle_contribution_at_R));
        R2 = imag(pz_non_only_imaginary) - real(pz_non_only_imaginary) * tan(-deg2rad(animationParams.min_angle_contribution_at_R));
        animationParams.R = max([abs(R1),abs(R2)]);
        %animationParams.R = 2*max([abs(animationParams.tf_zeros) abs(animationParams.tf_poles)]);     % radius of the D-curve. auto value is such that one can oftentimes see what happens in origin of nyquist curve
    end
    g = polesZeros2fncHandle(animationParams.tf_zeros, animationParams.tf_poles);
    
    %% preparations
    % prepare frame buffer
    animationParams.N_frames = animationParams.duration * animationParams.FPS;
    animationParams.N_trail = fix(animationParams.N_frames * animationParams.trail_length);
    
    % set up parametrisation according to duration and FPS
    animationParams.t_step = (animationParams.t_max-animationParams.t_min)/(animationParams.N_frames-1);
    t = animationParams.t_min:animationParams.t_step:animationParams.t_max;
    
    % prepare t-intervals that are plotted on each frame, considering the
    % set number of trail-values
    t_indexes = cell(animationParams.N_frames,1);
    for ii = 1:animationParams.N_frames
        t_indexes{ii} = max(1,ii-animationParams.N_trail):ii;
    end
    
    %% calculate the function values
    in = getInputFunctionHandle();
    %out = @(t,R,phi) g(in(t,animationParams));
    
    in_values = in(t,animationParams.R,animationParams.phi);
    [ts,in_values] = get_t_and_shape(animationParams);
    out_values = g(in_values);
    
    %% calculation of the plot(s)
    frames = drawFunctions(in_values, out_values, animationParams.tf_zeros, animationParams.tf_poles, animationParams.R, t, t_indexes, animationParams);
    
    if ~isstruct(frames) && frames == -1
        return;
    end
    
    %% display the movie player
    if ~animationParams.noPlayer && ~animationParams.export && exist('implay') == 2
        implay(frames,animationParams.FPS);
        set(findall(0,'tag','spcui_scope_framework'),'position',[150 150 length(frames(1).cdata(1,:,1))+20 length(frames(1).cdata(:,1,1))+30]);
    elseif ~(exist('implay') == 2)
        fprintf('Video player could not be found. Please install Matlabs'' ''Image Processing Toolbox'' in order to use this function\n');
    end
    
    % export the frames as a .mp4 file
    if animationParams.export
        if exist([animationParams.filename,'.mp4'],'file') == 2
            ii = 1;
            while exist([animationParams.filename,'_',num2str(ii),'.mp4'],'file') == 2
                ii = ii + 1;
            end
            filename = [animationParams.filename,'_',num2str(ii),'.mp4'];
        else
            filename = [animationParams.filename,'.mp4'];
        end
        v = VideoWriter(filename,'MPEG-4');
        v.FrameRate = animationParams.FPS;
        open(v);
        for ii = 1:length(frames)
            writeVideo(v,frames(ii));
        end
        close(v);
        fprintf(['Video successfully exported to "',filename,'"\n']);
    end
end

%% drawing function, returns a movie frame struct
function [frames] = drawFunctions(in_values, out_values, zeros, poles, R, t, t_indexes, animationParams)
    % remember how many zeros and poles we have
    zeros_N = length(zeros);
    poles_N = length(poles);

    % try to find optimal axis limits for the two plots, or let the user
    % decide if he chose to set either 'auto_size' to false
    [in_axis_xlim,in_axis_ylim,out_axis_xlim,out_axis_ylim] = autoZoom([zeros,poles],R,out_values);
    
    if ~animationParams.in_auto_size
       in_axis_xlim = [(animationParams.in_x0(1) - animationParams.in_dims(1)/2),(animationParams.in_x0(1) + animationParams.in_dims(1)/2)];
       in_axis_ylim = [(animationParams.in_x0(2) - animationParams.in_dims(2)/2),(animationParams.in_x0(2) + animationParams.in_dims(2)/2)];
    end
    
    if ~animationParams.out_auto_size
       out_axis_xlim = [(animationParams.out_x0(1) - animationParams.out_dims(1)/2),(animationParams.out_x0(1) + animationParams.out_dims(1)/2)];
       out_axis_ylim = [(animationParams.out_x0(2) - animationParams.out_dims(2)/2),(animationParams.out_x0(2) + animationParams.out_dims(2)/2)];
    end
    
    % derive the arrow length relative to subplot axis dimensions
    in_axis_width = in_axis_xlim(2) - in_axis_xlim(1);        % [1]
    in_axis_height = in_axis_ylim(2) - in_axis_ylim(1);       % [1]
    out_axis_width = out_axis_xlim(2) - out_axis_xlim(1);     % [1]
    out_axis_height = out_axis_ylim(2) - out_axis_ylim(1);    % [1]

    in_arrow_length = 0.04 * sqrt(in_axis_width^2 + in_axis_height^2);      % [1]
    out_arrow_length = 0.04 * sqrt(out_axis_width^2 + out_axis_height^2);   % [1]
    if out_arrow_length < 0.05
        out_arrow_length = 0.0001;
    end
    
    % preallocate space for frames to be stored
    frames(animationParams.N_frames) = struct('cdata',[],'colormap',[]);    
    
    try
        % create the figure with its two subplots and remember their handles
        fig = figure;
        sub1 = subplot(1,2,1);
        hold on;
        sub2 = subplot(1,2,2);
        hold on;

        % set the visual params of the subplots
        legacy = checkMatlabVersion();
        subplot(sub1);
        axis equal; % for 1:1 aspect ratio
        xlim(in_axis_xlim), ylim(in_axis_ylim);
        setaxes_origin(legacy);
        grid on;


        subplot(sub2);
        axis equal; % for 1:1 aspect ratio
        xlim(out_axis_xlim), ylim(out_axis_ylim);
        setaxes_origin(legacy);
        grid on;

        % calculate the correct positions (relative to interior of the figure)
        fig_plot_height = 2*animationParams.border + animationParams.plotSize;              % plot + border above and below, pixel
        fig_annotation_textbox_height = 14;                                                 % one box, pixel
        fig_annotation_height_sum = (max(zeros_N,poles_N)+2)*fig_annotation_textbox_height; % cumulative, pixel
        switch(legacy)
            case 'R2015b_or_newer'
                fig_legacy_correction = 0;
            otherwise
                fig_legacy_correction = 3*fig_annotation_textbox_height;
        end
        fig_height = fig_plot_height + fig_annotation_height_sum + fig_legacy_correction;   % pixel
        fig_width = 3*animationParams.border + 2*animationParams.plotSize;                                  % pixel

        % fig.Position expects pixels as unit
        fig.Position = [100 100 fig_width fig_height];

        % as subplot.Postion expects fractions of the inside of the figure, we
        % recalculate the pixel values
        fig_border_horizontal_frac = animationParams.border/fig_width;                      % fraction
        fig_border_vertical_frac = animationParams.border/fig_height;                       % fraction
        fig_plot_width_frac = (fig_width - 3*animationParams.border)/(2*fig_width);         % fraction
        fig_plot_height_frac = 1-(2*animationParams.border + fig_annotation_height_sum)/fig_height; % fraction

        % subX.Position expects fractions of the inside of the figure as unit
        sub1.Position = [fig_border_horizontal_frac, (1 - fig_border_vertical_frac - fig_plot_height_frac), fig_plot_width_frac, fig_plot_height_frac];
        sub2.Position = [(2*fig_border_horizontal_frac + fig_plot_width_frac),(1 - fig_border_vertical_frac - fig_plot_height_frac), fig_plot_width_frac, fig_plot_height_frac];

        % draw the zeros and poles in the left (input function) subplot
        subplot(sub1);
        for z = 1:zeros_N
            zero = zeros(z);
            zero_trunc = trunc(zero, in_axis_xlim, in_axis_ylim);
            if zero ~= zero_trunc
                % The zero in question being not equal to its truncated value
                % means that it lies outside the current subplot limits.
                % Instead of drawing just a circle at the border, plot an
                % arrow, indicating that there's an outlier.
                drawTextArrow([real(zero_trunc),imag(zero_trunc)],angle(zero),in_arrow_length,' o',[95 158 160]/255);
            else
                scatter(real(zero),imag(zero),60,'o','MarkerEdgeColor',[70 130 180]/255,'LineWidth',1.5);
            end
        end

        for p = 1:poles_N
            pole = poles(p);
            pole_trunc = trunc(pole, in_axis_xlim, in_axis_ylim);
            if pole ~= pole_trunc
                % same conclusion as above with the zeros
                drawTextArrow([real(pole_trunc),imag(pole_trunc)],angle(pole),in_arrow_length,' x',[255 140 0]/255);
            else
                scatter(real(pole),imag(pole),60,'x','MarkerEdgeColor',[255 140 0]/255,'LineWidth',1.5);
            end
        end

        % calculate the truncated values of the function (according to window
        % params)
        in_values_truncated = trunc(in_values, in_axis_xlim, in_axis_ylim);
        out_values_truncated = trunc(out_values, out_axis_xlim, out_axis_ylim);

        % plot the full input- and output curves with some transparency
        in_plot_full = plot(sub1,in_values_truncated,'Color',[0.05 0.4970 0.7410]); in_plot_full.Color(4) = 0.3;
        out_plot_full = plot(sub2,out_values_truncated,'Color',[0.05 0.4970 0.7410]); out_plot_full.Color(4) = 0.3;

        % prepare the curves' trail plots and remember their handle for later
        % use
        in_plot_trail = plot(sub1,0,0); set(in_plot_trail,'linewidth',2,'Color',[255 215 0]/255);
        out_plot_trail = plot(sub2,0,0); set(out_plot_trail,'linewidth',2,'Color',[255 215 0]/255);

        % prepare the arrows annotations that mark the value of the head of the
        % trail at the current frame and remember their handle for later use
        subplot(sub1); in_plot_arrow = annotation('line',[0 0],[1 0]);
        subplot(sub2); out_plot_arrow = annotation('line',[0 0],[1 0]);

        % We need to know the value of the head of the trail at the previous
        % frame in order to determine the arrow's direction. In the case of the
        % exact start direction, we need the value at the second to last frame
        % because if we took the last frame then we'd use f(t=1) = f(t=0).
        % For the debug case when t_min > 0 or t_max < 1 we simply
        % extrapolate from the second frame backwards.
        if t(1) == 0 && t(end) == 1
            in_values_head_prev = in_values(animationParams.N_frames - 1);
            out_values_head_prev = out_values(animationParams.N_frames - 1);
        else
            in_values_head_prev = 2*in_values(1) - in_values(2);
            out_values_head_prev = 2*out_values(1) - out_values(2);
        end

        % prepare the annotations below the plots and remember their handles

        % again we need to know relative (fractions) positions inside the figure
        fig_annotation_start_frac = (fig_annotation_height_sum + 10)/fig_height;    % fraction
        fig_annotation_textbox_frac = fig_annotation_textbox_height/fig_height;     % fraction

        figure(fig);
        % zero and pole contributions
        annotation('TextBox',[fig_border_horizontal_frac, (fig_annotation_start_frac - 0.5*fig_annotation_textbox_frac), fig_plot_width_frac, fig_annotation_textbox_frac],'String','Contribution of the zeros:','LineStyle','none','FontSize',9);
        for ii = 1:zeros_N
            z_annotations(ii) = annotation('TextBox',[fig_border_horizontal_frac, (fig_annotation_start_frac - fig_annotation_textbox_frac*(ii+1)), fig_plot_width_frac, fig_annotation_textbox_frac],'String',['zero ',num2str(ii)],'LineStyle','none','FontSize',9);
        end
        annotation('TextBox',[0.25, (fig_annotation_start_frac - 0.5*fig_annotation_textbox_frac), fig_plot_width_frac, fig_annotation_textbox_frac],'String','Contribution of the poles:','LineStyle','none','FontSize',9);
        for ii = 1:poles_N
            p_annotations(ii) = annotation('TextBox',[0.25, (fig_annotation_start_frac - fig_annotation_textbox_frac*(ii+1)), fig_plot_width_frac, fig_annotation_textbox_frac],'String',['pole ',num2str(ii)],'LineStyle','none','FontSize',9);
        end

        % cumulative calculations
        annotation('TextBox',[0.5 max(0,(fig_annotation_start_frac - 0.5*fig_annotation_textbox_frac)) fig_plot_width_frac fig_annotation_textbox_frac],'String','Resulting value of G:','LineStyle','none','FontSize',9);
        res_annotation1 = annotation('TextBox',[0.5 max(0,(fig_annotation_start_frac - 2*fig_annotation_textbox_frac)) fig_plot_width_frac fig_annotation_textbox_frac],'String','resultline 1','LineStyle','none','FontSize',9);
        res_annotation2 = annotation('TextBox',[0.5 max(0,(fig_annotation_start_frac - 3*fig_annotation_textbox_frac)) fig_plot_width_frac fig_annotation_textbox_frac],'String','resultline 2','LineStyle','none','FontSize',9);

        % main loop where we update the plots and annotations and store the
        % figure in a movie frame after each iteration
        for ii = 1:animationParams.N_frames        
            % update the trail plot data
            set(in_plot_trail,'XData',real(in_values_truncated(t_indexes{ii})),'YData',imag(in_values_truncated(t_indexes{ii})));
            set(out_plot_trail,'XData',real(out_values_truncated(t_indexes{ii})),'YData',imag(out_values_truncated(t_indexes{ii})));        

            % draw the (arrow-)head of the trail
            delete(in_plot_arrow); delete(out_plot_arrow);
            in_values_head = in_values(ii);
            out_values_head = out_values(ii);
            in_phi = angle(in_values_head - in_values_head_prev);
            out_phi = angle(out_values_head - out_values_head_prev);
            subplot(sub1); in_plot_arrow = drawArrow([real(in_values_truncated(ii)),imag(in_values_truncated(ii))],in_phi,in_arrow_length);
            subplot(sub2); out_plot_arrow = drawArrow([real(out_values_truncated(ii)),imag(out_values_truncated(ii))],out_phi,out_arrow_length);

            % update the title
            title(sub1,['Current value of D-Curve: ',num2str(in_values_head,'%.1f'),': M = ',num2str(abs(in_values_head),'%.2f'),' p = ',num2str(rad2deg(angle(in_values_head)),'%.2f'),'°']);
            title(sub2,['Nyquist: G(',num2str(in_values_head,'%.1f'),') = ',num2str(out_values_head,'%.1f'),': M = ',num2str(abs(out_values_head),'%.2f'),' p = ',num2str(rad2deg(angle(out_values_head)),'%.2f'),'°']);

            % update the zero- and pole contribution and the cumulative values
            res_magnitude = '(';
            res_phase = '';
            for z = 1:zeros_N
                z_contribution = in_values(ii) - zeros(z);
                z_annotations(z).String = ['Z',num2str(z),': (',num2str(in_values(ii),'%.1f'),') - (',num2str(zeros(z),'%.1f'),'): M=',num2str(abs(z_contribution),'%.1f'),' p=',num2str(rad2deg(angle(z_contribution)),'%.1f'),'°'];
                res_magnitude = [res_magnitude,'*',num2str(abs(z_contribution),'%.2f')];
                res_phase = [res_phase,'+',num2str(rad2deg(angle(z_contribution)),'%.2f')];
            end
            res_magnitude = [res_magnitude,')/('];
            for p = 1:poles_N
                p_contribution = in_values(ii) - poles(p);
                p_annotations(p).String = ['P',num2str(p),': (',num2str(in_values(ii),'%.1f'),') - (',num2str(poles(p),'%.1f'),'): M=',num2str(abs(p_contribution),'%.1f'),' p=',num2str(rad2deg(angle(p_contribution)),'%.1f'),'°'];            
                res_magnitude = [res_magnitude,'*',num2str(abs(p_contribution),'%.2f')];
                res_phase = [res_phase,'-',num2str(rad2deg(angle(p_contribution)),'%.2f')];
            end
            res_magnitude = [res_magnitude,') = ',num2str(abs(out_values(ii)),'%.3f')];
            res_phase = [res_phase,') = ',num2str(rad2deg(angle(out_values(ii))),'%.3f')];

            res_annotation1.String = ['Magnitude: ',res_magnitude];
            res_annotation2.String = ['Phase:       ',res_phase];

            % save the current frame
            frames(ii) = getframe(gcf);

            in_values_head_prev = in_values_head;
            out_values_head_prev = out_values_head;       
            
            %keyboard;
        end
    catch err
        if strcmp(err.identifier,'MATLAB:structRefFromNonStruct')
            frames = -1;
            return;
        else
            rethrow(err);
        end
    end
    if exist('implay','file') == 2
        close(fig);
    end
end

% The purpose of this function is to provide a means of detecting clumping
% of zeros and poles. If we have a bunch of zeros and poles that are fairly
% close to each other with one outlier, forget about the outlier and zoom
% in on the group.
function [in_axis_xlim,in_axis_ylim,out_axis_xlim,out_axis_ylim] = autoZoom(zp, R, out_values)
    % find the 'center of gravity' of the zeros and poles
    locus_m = mean(zp);
    
    % where are 66% of the points?
    locus_std_real = std(real(zp));
    locus_std_imag = std(imag(zp));
    
    % get indexes of those points that lie within the ellipse given by the stddev in x- and y-direction
    locus_inner =  find(abs(zp-locus_m) <= sqrt(locus_std_real^2+locus_std_imag^2));
    
    % how much area do these points cover in the zoomed-out plot (when the D-curve fits into it)?
    locus_inner_dims = [max(real(zp(locus_inner)))-min(real(zp(locus_inner)));max(imag(zp(locus_inner)))-min(imag(zp(locus_inner)))];
    if min(locus_inner_dims)/max(locus_inner_dims) <= 1/50
        % almost one-dimensional
        locus_inner_dim = max(locus_inner_dims);
        locus_fract = locus_inner_dim/(2.01*R);
    else
        locus_inner_area = locus_inner_dims(1)*locus_inner_dims(2);
        locus_fract = locus_inner_area/(2.01*R)^2;
    end
    
    % if that area is below 0.2% then we zoom into the bulk of the points,
    % leaving out outliers
    if locus_fract > 0.015
        % all good, fit the whole D-curve and the zp-points into the plot
        in_axis_xlim = 1.05*R*[-1,1];
        in_axis_ylim = 1.05*R*[-1,1];
    else
        % include the origin in the box around the clumped points and then
        % make the plot size 8 times (rule of thumb) larger
        in_axis_xlim = stretch_centered([min([real(zp(locus_inner)),0]),max([real(zp(locus_inner)),0])],8);
        in_axis_ylim = stretch_centered([min([imag(zp(locus_inner)),0]),max([imag(zp(locus_inner)),0])],8);
        
        width = in_axis_xlim(2) - in_axis_xlim(1);
        height = in_axis_ylim(2) - in_axis_ylim(1);
        
        % make the plot rectangular (make the smaller axis as large as the
        % bigger one)
        if width >= height
            y0 = mean(in_axis_ylim);
            in_axis_ylim = [(y0-width/2),(y0+width/2)];
            height = width;
        else
            x0 = mean(in_axis_xlim);
            in_axis_xlim = [(x0-height/2),(x0+height/2)];
            width = height;
        end
    end
    
    % for the right subplot it's fairly easy. Have the Plot such that the
    % whole nyquist curve fits into it.
    out_axis_xlim = stretch_centered([min([real(out_values);0]),max([real(out_values);0])],1.15);
    out_axis_ylim = stretch_centered([min([imag(out_values);0]),max([imag(out_values);0])],1.15);
    
    width = out_axis_xlim(2) - out_axis_xlim(1);
    height = out_axis_ylim(2) - out_axis_ylim(1);

    % make this plot also rectangular
    if width >= height
        y0 = mean(out_axis_ylim);
        out_axis_ylim = [(y0-width/2),(y0+width/2)];
        height = width;
    else
        x0 = mean(out_axis_xlim);
        out_axis_xlim = [(x0-height/2),(x0+height/2)];
        width = height;
    end

end

function arrowHandle = drawArrow(x0, phi, l)
    r = l*[cos(phi),sin(phi)];
    
    arrowHandle = annotation('Arrow');
    set(arrowHandle,'parent',gca,'position',[x0-r,r]);
end

function arrowHandle = drawTextArrow(x0, phi, l, text, color)
    r = l*[cos(phi),sin(phi)];
    
    arrowHandle = annotation('TextArrow');
    set(arrowHandle,'parent',gca,'position',[x0-r,r],'String',text,'Color',color);
end

%% input function shapes
% t is in the following range 0<=t<=1.
% the shape should be continuous and a closed curve (i.e. shape(0)=shape(1))
function inputFunctionHandle = getInputFunctionHandle()
    %inputFunctionHandle = @(t) unitCircle(t);
    %inputFunctionHandle = @(t) dCurve(t);
    inputFunctionHandle = @(t,R,phi) roundedCurve(t,R,phi);
end

% this D-shaped curve has small arcs at the edges to smooth them and
% make the derivative of the curve continuous
% The shape is basically split into 5 parts to which the indexes correspond
% The two parameters one can play with are phi, which determines the degree
% of roundoff (0<phi<180) and the radius of the big arc R.
function [t, shape] = get_t_and_shape(animationParams)
    poles_zeros = [animationParams.tf_zeros animationParams.tf_poles];
    phi = animationParams.phi;
    R = animationParams.R;
    
    detour_interval_max = 1/2;
    detour_angle = 45;
    
    t_min = animationParams.t_min;
    t_max = animationParams.t_max;
    
    assert(t_min >= 0 && t_min <1 && t_max > 0 && t_max <= 1 && t_min < t_max);
    t_step = (t_max-t_min)/(animationParams.N_frames-1);
    t = t_min:t_step:t_max;
    
    n_steps = length(t);
    
    phi = deg2rad(phi);
    
    % calculate values of the small round-off circles
    r = R*sin(phi)/(1+sin(phi));
    m2 = r + 1i*sqrt(R^2-2*R*r);
    m4 = r - 1i*sqrt(R^2-2*R*r);
    
    % need to know how long the parts are
    d1 = sqrt(R^2-2*R*r);
    d2 = r*(pi/2 + phi);
    d3 = R*(pi - 2*phi);
    d4 = d2;
    d5 = d1;
    
    % length of the round part (small arcs + big arc)
    d_arcs = d2 + d3 + d4;
    
    % calculate the t-values at which two adjacent regions touch
    % the only arbitrary choice here is t12. together with t45 it
    % determines the fraction of time we need for the y-axis part of the
    % D-curve.
    % the v-values are velocity-factors
    t12 = 2/5; 
    t45 = 1-t12;
    v_ave = d_arcs / (t45 - t12);
    
    t23 = t12 + d2/v_ave;
    t34 = 1-t23;
    
    v11 = 3; % again the only arbitrary choice amongst the velocity-factors. determines the speed of the animation near the origin.
    v12 = (d1-v11*t12)/t12^4;
    v51 = v11;
    v52 = v12;
        
    ts = t_min:t_step:t_max;    % '1st row': t-values
    shape = zeros(n_steps,1);   % '2nd row': shape-values
    
    region1_idx = find(ts >=t_min          & ts <= min(t_max,t12));
    region2_idx = find(ts > max(t_min,t12) & ts <= min(t_max,t23));
    region3_idx = find(ts > max(t_min,t23) & ts <= min(t_max,t34));
    region4_idx = find(ts > max(t_min,t34) & ts <= min(t_max,t45));
    region5_idx = find(ts > max(t_min,t45) & ts <=t_max);
    
    shape(region1_idx) = (v11*ts(region1_idx) + v12*ts(region1_idx).^4)*1i; % part 1
    shape(region2_idx) = m2 + r*exp(1i*(pi - (ts(region2_idx)-t12).*(pi/2 + phi)/(t23-t12))); % part 2
    shape(region3_idx) = R*exp(1i*(pi/2 - phi - (ts(region3_idx)-t23).*(pi-2*phi)/(t34-t23))); % part 3
    shape(region4_idx) = m4 + r*exp(1i*(phi - pi/2 - (ts(region4_idx)-t34).*(pi/2 + phi)/(t45-t34))); % part 4
    shape(region5_idx) = -(v51*(1-ts(region5_idx)) + v52*(1-ts(region5_idx)).^4)*1i; % part 5
    
    % get poles and zeros which lie on the imaginary axis
    imag_pz_idx = find(real(poles_zeros) == 0);
    if isempty(imag_pz_idx)
        return;
    end
    
    imaginary_pz = poles_zeros(imag_pz_idx);
    
    % find minimum distance between any pair of only-imaginary poles and
    % zeros
    mask = triu(ones(length(imaginary_pz)),1);
    [pz1,pz2] = meshgrid(imaginary_pz,imaginary_pz);
    pairwise_distances = abs(pz1.*mask - pz2.*mask);
    d_min = min(nonzeros(pairwise_distances));
    if isempty(d_min)
        d_min = 1e200;
    end
    
    % determine the radii for the detours 
    % 1)    based on the calculated minimum distances between imaginary
    %       poles and zeros
    % 2)    based on the provided maximal secant length
    detour_angle = detour_angle*pi/180;
    max_pairwise_detour_radius = sqrt((d_min/2)^2/(1-tan(detour_angle)^2/(4+tan(detour_angle)^2)));
    max_parametrized_detour_radius =  sqrt((detour_interval_max/2)^2/(1-tan(detour_angle)^2/(4+tan(detour_angle)^2)));
    
    % set the parameters for all of the detour circles
    secant_detour = min(d_min,detour_interval_max);
    radius_detour = min(max_pairwise_detour_radius,max_parametrized_detour_radius);
    x_shift_detour = sqrt(radius_detour^2 - (secant_detour/2)^2);
    
    for ii = 1:length(imaginary_pz)
        detour_interval_low = imag(imaginary_pz(ii))-secant_detour/2;
        detour_interval_high = imag(imaginary_pz(ii))+secant_detour/2;
        detour_point_idx = find(real(shape) == 0 & imag(shape) > detour_interval_low & imag(shape) < detour_interval_high);
        shape(detour_point_idx) = shifted_circle(shape(detour_point_idx),imaginary_pz(ii));
    end
    
    
    function z = shifted_circle(t,y_shift_detour)
        z = sqrt(radius_detour^2-imag(t-y_shift_detour).^2) - x_shift_detour + imag(t).*1i;
    end
end

%% transfer function preparation
function g = polesZeros2fncHandle(zeros,poles)
    g = @(z) polyval(poly(zeros),z)./polyval(poly(poles),z);
end

%% various stuff
% this function clips the input values to the maximum value allowed inside
% the plot
function z = trunc(x,xlim,ylim)
    z = max(xlim(1), min(xlim(2), real(x))) ...
        + 1i*max(ylim(1), min(ylim(2), imag(x)));
end

% stretches a given interval with a given factor about its center
function result = stretch_centered(in, factor)
    size_in = in(2) - in(1);
    result = [(in(1)-(factor-1)*size_in/2),(in(2)+(factor-1)*size_in/2)];
end

function legacy = checkMatlabVersion()
    % release history and numbers:
    % https://en.wikipedia.org/wiki/MATLAB#Release_history
    
	% https://ch.mathworks.com/help/doc-archives.html
	
    % whole release download:
    % https://ch.mathworks.com/downloads/web_downloads/select_release
    
    % only runtime download:
    % https://ch.mathworks.com/products/compiler/mcr/
    
    v = ver;
    idx = find(ismember({v(:).Name},'MATLAB'));
    matlab_version = str2double(v(idx).Version(1:3));

    if matlab_version <= 8.3
        % we want at least HG2 to work with.
        legacy = 'too_old';
        if ~exist('animated_nyquist_plot_warningOff','file')
            fprintf('Warning: You are using a release of Matlab that doesn''t support the current\n')
            fprintf('HG2 graphics engine. Some features such as axle repositioning are deactivated, but there''s still a chance that the script crashes.\n')
            fprintf('Please consider upgrading to R2014b or higher.\n');
            fprintf('This script has not been tested with releases lower than R2014b.\n');
            fprintf('(Suppress this warning by placing an empty file named ''animated_nyquist_plot_warningOff'' in the script''s working folder.)\n');
        end
    elseif matlab_version < 8.6
        legacy = 'pre_R2015b';
    else
        legacy = 'R2015b_or_newer';
    end
    
    if ~exist('animated_nyquist_plot_warningOff','file') && matlab_version >= 8.5 && matlab_version <= 9.0
        fprintf('Note: this script uses release-specific functions of Matlab R2014b\n');
        fprintf('respectively from R2015b on. Testing has been done only for R2014b and R2016b though.\n');
        fprintf('Please report any bugs together with their input!\n');
        fprintf('(Suppress this warning by placing an empty file named ''animated_nyquist_plot_warningOff'' in the script''s working folder.)\n');
    end
end

% moves the current ax(l)es into the origin of the plot
function [] = setaxes_origin(legacy)
    switch(legacy)
        case 'R2015b_or_newer'
            ax = gca;
            ax.XAxisLocation = 'origin';
            ax.YAxisLocation = 'origin';
        case 'pre_R2015b'
            % http://undocumentedmatlab.com/blog/customizing-axes-part-2
            ax = gca;
            ax.XBaseline.Color = 'k';
            ax.YBaseline.Color = 'k';
            
            ax.XBaseline.LineWidth = 1.5;
            ax.YBaseline.LineWidth = 1.5;
            
            ax.YBaseline.Visible = 'on';
            ax.XBaseline.Visible = 'on';
            
            ax.XRuler.Axle.Visible = 'off';
            ax.YRuler.Axle.Visible = 'off';

            ax.LineWidth = 1;

            % 2nd possibility for workaround: http://matlabvn.com/matlab-plotting-tricks-move-x-axis-center/
            % other:
            % https://ch.mathworks.com/matlabcentral/answers/147087-how-to-change-the-axes-position-in-matlab
            % http://ch.mathworks.com/matlabcentral/fileexchange/22956-axescenter
        otherwise
            % Too old to produce compatible code. Just skip the axle
            % repositioning.
            return;
    end
end