% -------------------------------------------------------------------------
%   ANP: Animated Nyquist Plot
%   **************************
%
%   Author: Stefan Rickli, stefanrickli [at] gmx.ch
%           https://blogs.ethz.ch/ricklis
%           https://github.com/StefanRickli/anp
%
%   Version 5.3.0
%   
%   
%   What does it do?:
%   -----------------
%   This script calculates for a given set of zeros and poles a suitable
%   D-curve, plots it together with the zeros and poles and draws the
%   corresponding nyquist curve in another plot.
%   For every time step one is shown the contribution of each zero and pole
%   to the magnitude and phase of the nyquist curve's values.
%   One can also go forward and backward in the animation by using the
%   appropriate buttons in the toolbar.
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
%           D-curve manually. Per default it choses its size such that the
%           every pole or zero at least contributes 65° in the nyquist plot
%           when the straight part of the D-curve goes into the
%           'infinity'-half circle.
%   
%   'plot_size':    Specifies the dimensions of each of the plots in
%                   pixels.
%                   Default size is 500px.
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
%   ANP can't handle input- or output delays and other nonlinearities, yet.
%
%
%   Known Issues:
%   -------------
%
%   Please visit https://github.com/StefanRickli/anp/issues to find a list
%   of documented issues.
%
%
%   Announcements:
%   --------------
%   The next major version will see basic MIMO support.
%   Along with it the D-contour will see a change, i.e. detours around
%   poles/zeros that lie exactly on the Im-axis will go around them on the
%   LHP and not the RHP anymore. This has no impact for the Nyquist
%   stability criterion in the SISO case, but is very vital to the MIMO
%   Nyquist stability criterion.
%
% -------------------------------------------------------------------------

function arg_out = anp_main(varargin)
    addpath('anp_files');
    addpath('anp_files/mimo');
    addpath('anp_files/d_contour_tiny_fcts');
    addpath('anp_files/tools');
    addpath('anp_icons');
    
    anp_check_Matlab_version(); % sets global variable 'matlab_version'
    
    global debug_graphics debug_graphics_interpolation debug_text;
    debug_graphics = false;
    debug_graphics_interpolation = false;
    debug_text = false;
    debug_no_reuse = false;
    
    args = anp_parse_arguments(varargin{:});
    
    % if argument checking has passed, get previous GUI and TF-Processing
    % objects back or re-initialize them if they have been deleted in the
    % meantime
    persistent h_tf_processor
    
    try
        h_gui = anp_gui;
        
        if isempty(h_tf_processor) || ~isvalid(h_tf_processor) || debug_no_reuse
            h_tf_processor =        tf_processor;
        else
            fprintf('(Reusing old tf_processor instance. If subsequent calls to anp_main behave odd, try ''clear all'')\n\n');
        end
        
        % pass the handles to the GUI and tf_processor in the args
        args.gui_handle =           h_gui;
        args.processor_handle =     h_tf_processor;
        
        args.tf_poles =             roots(args.tf_obj.Denominator{1}).';
        args.tf_zeros =             roots(args.tf_obj.Numerator{1}).';
        
        % this kicks off the calculation of the data that the
        % GUI object then displays
        h_tf_processor.init_params(args);
        
        % let the GUI read the relevant parameters from the arguments
        h_gui.init_params(args);
        
        h_gui.init_visuals();
        
        % Handling of debugging keywords
        if args.trigger_step
            h_gui.trigger_step();
        end
        
        if args.return_handle
            arg_out = {h_tf_processor, h_gui};
        end
        
    catch err
        clear h_tf_processor;
        disp(getReport(err,'extended'));
        
        if args.cleanup_after_error
            delete(h_gui);
        end
    end
end

