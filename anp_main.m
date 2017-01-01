% -------------------------------------------------------------------------
%   Animated Nyquist Plot
%   *********************
%
%   Author: Stefan Rickli, ricklis [at] student.ethz.ch
%           https://blogs.ethz.ch/ricklis
%
%   Version 5.0, 2017-01-01
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
%   'duration': deprecated, will be ignored.
%   
%   'FPS':      deprecated, will be ignored.
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
%   1)  The function can't handle input- or output delays.
%   2)  v5.0 has only been minimally tested with optional arguments.
%
%
%   Known Issues:
%   -------------
%
%   - Pole and zero contributions of those with positive real part may be
%     incorrect.
%   - Transfer functions with a high relative degree of numerator and
%     denominator produce plots with very small amplitudes which 
%     lead Matlab to draw arrow-annotations incorrectly.
%
% -------------------------------------------------------------------------

function [] = anp_main(varargin)
    
    global debug_graphics debug_text;
    debug_graphics = false;
    debug_text = false;

    addpath('anp_files');
    addpath('anp_icons');

    args =                      anp_parse_arguments(varargin{:});

    % if argument checking has passed, get previous GUI and TF-Processing
    % objects back or re-initialize them if they have been deleted in the
    % meantime
    persistent h_tf_processor

    h_gui = anp_gui;

    if isempty(h_tf_processor) || ~isvalid(h_tf_processor)
        h_tf_processor =        anp_tf_processor;
    else
        fprintf('Reusing old TF processor instance\n');
    end

    args.gui_handle =           h_gui;
    args.processor_handle =     h_tf_processor;
    args.tf_poles =             roots(args.tf_obj.Numerator{1})';
    args.tf_zeros =             roots(args.tf_obj.Denominator{1})';

    h_tf_processor.init_params(args);

    h_gui.init_params(args);

    h_gui.init_visuals();

end

